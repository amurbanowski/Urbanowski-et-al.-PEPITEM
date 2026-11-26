#monocyte and macrophage subclustering, and a final population/exclusion reason for every cell

```{r}

library(Seurat)
library(dplyr)

dir <- "/rds/projects/m/mcgetthm-rnaseq-arthritis/AU_final_analysis"
inp <- file.path(dir, "AU_input_objects")
out <- file.path(dir, "AU_final_figures_update", "exclusion_flow_sensitivity")
dir.create(out, recursive = T, showWarnings = F)

aggr <- readRDS(file.path(inp, "all_synovial_object_used.rds"))
DefaultAssay(aggr) <- "RNA"

#same steps for every subclustering, on the RNA assay
recluster <- function(obj, npcs = 30, res) {
  if (inherits(obj[["RNA"]], "Assay5")) obj <- JoinLayers(obj)
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, nfeatures = 2000)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, npcs = npcs)
  obj <- FindNeighbors(obj, dims = 1:20)
  obj <- FindClusters(obj, resolution = res)
  obj <- RunUMAP(obj, dims = 1:20)
  obj
}

```

```{r}

#monocytes, 9529 cells
mono <- subset(aggr, Cluster_all == "monocytes")
mono <- recluster(mono, res = c(0.1, 0.2, 0.3, 0.5))
mono$Monocyte_subcluster <- mono$RNA_snn_res.0.3

DimPlot(mono, group.by = "Monocyte_subcluster", label = T)
saveRDS(mono, file.path(inp, "monocyte_subclustered_object.rds"))

#cluster 0 re-clustered individually to pull out the stress cells, 3520 cells
c0 <- subset(mono, Monocyte_subcluster == "0")
c0 <- recluster(c0, npcs = 20, res = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1))
c0$Cluster0_subcluster <- c0$RNA_snn_res.0.5

DimPlot(c0, group.by = "Cluster0_subcluster", label = T)

```

```{r}

#macrophages first pass, 11059 cells
mac <- subset(aggr, Cluster_all %in% c("M1 macrophages", "M2 macrophages"))
mac <- recluster(mac, res = c(0.2, 0.3, 0.4, 0.6, 0.8))
mac$Macrophage_subcluster <- mac$RNA_snn_res.0.4

DimPlot(mac, group.by = "Macrophage_subcluster", label = T)
saveRDS(mac, file.path(inp, "macrophage_subclustered_object.rds"))

#cleaned macrophages, first pass 0/1/4/6/8/9, 6331 cells
mac_clean <- subset(mac, Macrophage_subcluster %in% c("0", "1", "4", "6", "8", "9"))
mac_clean <- recluster(mac_clean, res = c(0.1, 0.15, 0.2, 0.3, 0.4))
mac_clean$Macrophage_subcluster_cleaned <- mac_clean$RNA_snn_res.0.2

DimPlot(mac_clean, group.by = "Macrophage_subcluster_cleaned", label = T)
saveRDS(mac_clean, file.path(inp, "macrophage_subclustered_object_cleaned.rds"))

#first pass 2/3 (neutrophil genes) re-clustered to separate neutrophils from macrophages, 2921 cells
neut <- subset(mac, Macrophage_subcluster %in% c("2", "3"))
neut <- recluster(neut, res = 0.4)
neut$Neutrophil_recluster_subcluster <- neut$RNA_snn_res.0.4

DimPlot(neut, group.by = "Neutrophil_recluster_subcluster", label = T)

```

#labels (from the marker genes for each subcluster)

```{r}

#monocyte subclusters
mono_pop <- c("0" = "Classical", "5" = "Transitional", "6" = "Non-classical")
mono_why <- c("1" = "Neutrophil contamination (mature granule genes)",
              "2" = "Neutrophil-lineage contamination",
              "3" = "Mitochondrial / dissociation stress",
              "4" = "Ambiguous low-confidence myeloid")

#cluster 0: subcluster 0 is Atf3/Jun/Ccl3/Ccl4/Hspa4l stress, 1-5 are classical
c0_stress <- colnames(c0)[c0$Cluster0_subcluster == "0"]

#cleaned macrophage subclusters
mac_pop <- c("0" = "MHC-II-high", "1" = "Classical", "2" = "TLF+ resident-like", "5" = "TLF+ resident-like",
             "6" = "TLF+ resident-like")
mac_why <- c("3" = "Stress / dissociation (macrophage subclustering)", "4" = "Neutrophils (macrophage subclustering)")

#first pass subclusters that were not macrophages
mac_moved <- c("10" = "Non-classical", "5" = "Dendritic cells", "7" = "Dendritic cells", "12" = "Dendritic cells",
               "13" = "Dendritic cells", "11" = "B cells")

#neutrophil re-cluster: 2/4/5/6/7 are neutrophils, 0/1/3 looked like macrophages but were not carried forward
neut_cells <- colnames(neut)[neut$Neutrophil_recluster_subcluster %in% c("2", "4", "5", "6", "7")]

```

#final labels

```{r}

fate <- data.frame(Barcode = colnames(aggr), Mouse = aggr$Sample_all, Condition = aggr$Condition_all,
                   Cluster_all = aggr$Cluster_all)

sample_to_condition <- c(S1 = "Control", S5 = "Control", S3 = "PEPITEM", S4 = "PEPITEM", S7 = "PEPITEM", S8 = "PEPITEM")
stopifnot(all(fate$Condition == sample_to_condition[fate$Mouse]))
fate$Compartment <- "Other compartments"
fate$Compartment[fate$Cluster_all == "monocytes"] <- "Monocyte compartment"
fate$Compartment[fate$Cluster_all %in% c("M1 macrophages", "M2 macrophages")] <- "Macrophage compartment"
fate$Final_population <- NA
fate$Exclusion_reason <- NA
fate$Reassigned_between_compartments <- FALSE

#other cell types are kept as they are
k <- fate$Compartment == "Other compartments"
x <- sub("myloid", "myeloid", fate$Cluster_all[k])
fate$Final_population[k] <- paste0(toupper(substr(x, 1, 1)), substring(x, 2))

#monocytes
i <- which(fate$Compartment == "Monocyte compartment")
s <- as.character(mono$Monocyte_subcluster[match(fate$Barcode[i], colnames(mono))])
fate$Final_population[i] <- unname(mono_pop[s])
fate$Exclusion_reason[i] <- unname(mono_why[s])
k <- i[fate$Barcode[i] %in% c0_stress]
fate$Final_population[k] <- NA
fate$Exclusion_reason[k] <- "Dissociation-stress artifact (monocyte subcluster 0)"

#macrophages
j <- which(fate$Compartment == "Macrophage compartment")
bc <- fate$Barcode[j]
orig <- as.character(mac$Macrophage_subcluster[match(bc, colnames(mac))])
cl <- as.character(mac_clean$Macrophage_subcluster_cleaned[match(bc, colnames(mac_clean))])

pop <- unname(mac_pop[cl])
why <- unname(mac_why[cl])
moved <- !is.na(cl) & cl == "1"

k <- orig %in% names(mac_moved)
pop[k] <- unname(mac_moved[orig[k]])
moved[k] <- TRUE

k <- orig %in% c("2", "3")
why[k] <- ifelse(bc[k] %in% neut_cells, "Neutrophils (confirmed by reclustering subclusters 2/3)",
                 "Subclusters 2/3 flagged retain-macrophage but not carried forward")

fate$Final_population[j] <- pop
fate$Exclusion_reason[j] <- why
fate$Reassigned_between_compartments[j] <- moved

fate$Status <- ifelse(is.na(fate$Exclusion_reason), "Retained", "Excluded")
write.csv(fate, file.path(out, "cell_fate_per_cell.csv"), row.names = F)

#should be 22108 retained and 10845 excluded, with Classical 3702, Transitional 191, Non-classical 263,
#MHC-II-high 2214, TLF+ resident-like 1713, Dendritic cells 1515 (all moved from the macrophage compartment)
table(fate$Status)
table(fate$Final_population[fate$Status == "Retained"])

#exclusions per mouse
table(fate$Exclusion_reason, fate$Mouse)

```
