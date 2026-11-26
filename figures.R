#UMAPs, marker dot plots and functional gene panels
#needs cell_fate_per_cell.csv from subclustering.R

```{r}

library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)

dir <- "/rds/projects/m/mcgetthm-rnaseq-arthritis/AU_final_analysis"
inp <- file.path(dir, "AU_input_objects")
out <- file.path(dir, "AU_figures_update")
dir.create(out, recursive = T, showWarnings = F)

fate <- read.csv(file.path(dir, "AU_final_figures_update", "exclusion_flow_sensitivity",
                           "cell_fate_per_cell.csv"))
fate <- fate[fate$Status == "Retained", ]

aggr <- readRDS(file.path(inp, "all_synovial_object_used.rds"))
DefaultAssay(aggr) <- "RNA"
aggr <- aggr[, fate$Barcode]
if (inherits(aggr[["RNA"]], "Assay5")) aggr <- JoinLayers(aggr)
aggr$population <- fate$Final_population[match(colnames(aggr), fate$Barcode)]

#broad populations - monocyte and macrophage subsets collapsed back together
mono_types <- c("Classical", "Transitional", "Non-classical")
mac_types <- c("MHC-II-high", "TLF+ resident-like")
aggr$broad <- aggr$population
aggr$broad[aggr$population %in% mono_types] <- "Monocytes"
aggr$broad[aggr$population %in% mac_types] <- "Macrophages"

table(aggr$broad)

```

```{r}

#UMAPs - whole dataset by broad population, subsets on their own subclustering UMAPs
p <- DimPlot(aggr, group.by = "broad") + NoAxes()
ggsave(file.path(out, "umap_broad.pdf"), p, width = 6, height = 4.5)

mono_umap <- readRDS(file.path(inp, "monocyte_subclustered_object.rds"))
mono_umap <- mono_umap[, intersect(colnames(mono_umap), colnames(aggr)[aggr$population %in% mono_types])]
mono_umap$population <- factor(aggr$population[match(colnames(mono_umap), colnames(aggr))], levels = mono_types)
p <- DimPlot(mono_umap, group.by = "population") + NoAxes()
ggsave(file.path(out, "umap_monocytes.pdf"), p, width = 5, height = 4)

mac_umap <- readRDS(file.path(inp, "macrophage_subclustered_object_cleaned.rds"))
mac_umap <- mac_umap[, intersect(colnames(mac_umap), colnames(aggr)[aggr$population %in% mac_types])]
mac_umap$population <- factor(aggr$population[match(colnames(mac_umap), colnames(aggr))], levels = mac_types)
p <- DimPlot(mac_umap, group.by = "population") + NoAxes()
ggsave(file.path(out, "umap_macrophages.pdf"), p, width = 5, height = 4)
rm(mono_umap, mac_umap)

```

```{r}

#marker dot plots
broad_markers <- c("Adgre1", "C1qa", "C1qb", "C1qc", "Fcgr1", "Mertk", "Timd4", "Lyve1", "Folr2", "Vsig4",
                   "Ly6c2", "Ccr2", "Sell", "Plac8", "Cx3cr1", "Nr4a1",
                   "Flt3", "Itgax", "H2-Ab1", "Cd74", "Xcr1", "Clec10a",
                   "Cd3d", "Cd3e", "Trbc1", "Cd4", "Cd8a",
                   "Cd79a", "Cd79b", "Ms4a1", "Cd19", "Cd22",
                   "Cpa3", "Kit", "Tpsb2", "Fcer1a", "Mcpt4", "Mcpt8",
                   "Mki67", "Top2a", "Tubb5", "Stmn1", "Birc5",
                   "Col1a1", "Dcn", "Pdgfra", "Pecam1", "Vwf", "Hba-a1")
mono_markers <- c("Ly6c2", "Ccr2", "Sell", "Chil3", "Plac8", "Vcan", "Cx3cr1", "Nr4a1", "Ace", "Itgal",
                  "Adgre4", "Spn", "Csf1r", "Lyz2", "Cd14", "Fcgr1")
mac_markers <- c("Cx3cr1", "Vsig4", "Folr2", "Timd4", "Lyve1", "H2-Aa", "H2-Ab1", "H2-Eb1", "Cd74", "Ciita",
                 "Il1b", "Tnf", "Ccr2", "S100a8", "S100a9", "Cxcl9", "Cxcl10", "Spp1", "Mmp12", "Mmp14", "Fn1",
                 "Trem2", "Mrc1", "Mertk", "Axl", "Gas6", "Stab1", "Csf1r", "C1qa", "C1qb", "C1qc", "Adgre1",
                 "Fcgr1", "Mafb")

dot_colours <- scale_colour_gradient2(low = "#0057D9", mid = "#FFFFFF", high = "#E31A1C", midpoint = 0)

p <- DotPlot(aggr, features = broad_markers, group.by = "broad") + dot_colours + RotatedAxis()
ggsave(file.path(out, "dotplot_broad.pdf"), p, width = 14, height = 4)

mono <- aggr[, aggr$population %in% mono_types]
mono$population <- factor(mono$population, levels = mono_types)
p <- DotPlot(mono, features = mono_markers, group.by = "population") + dot_colours + RotatedAxis()
ggsave(file.path(out, "dotplot_monocytes.pdf"), p, width = 7, height = 2.5)

mac <- aggr[, aggr$population %in% mac_types]
mac$population <- factor(mac$population, levels = mac_types)
p <- DotPlot(mac, features = mac_markers, group.by = "population") + dot_colours + RotatedAxis()
ggsave(file.path(out, "dotplot_macrophages.pdf"), p, width = 11, height = 2.2)

```

```{r}

#functional gene panels (Table S6 and S7)
mono_genes <- list(
  "Inflammatory cytokines / alarmins" = c("Il1b", "Tnf", "Il6", "Ptgs2", "S100a8", "S100a9"),
  "Inflammasome / pyroptosis" = c("Nlrp3", "Pycard", "Casp1", "Gsdmd"),
  "Interferon response" = c("Stat1", "Irf1", "Isg15", "Ifit1", "Ifit3", "Gbp2"),
  "Chemokine secretion / recruitment" = c("Ccl2", "Ccl3", "Ccl4", "Ccl7", "Cxcl2", "Cxcl9", "Cxcl10"),
  "Chemokine sensing / trafficking" = c("Ccr1", "Ccr2", "Ccr5", "Cxcr4"),
  "Adhesion / transendothelial migration" = c("Sell", "Selplg", "Itgam", "Itgb2", "Cd44", "Lsp1", "Rac2"),
  "Patrolling / vascular surveillance" = c("Nr4a1", "Cx3cr1", "Itgal", "Treml4", "Cd300a", "Lair1"),
  "Antigen processing / presentation" = c("H2-Aa", "H2-Ab1", "H2-Eb1", "Cd74", "Ciita", "Ctss"),
  "Phagocytosis / lysosome" = c("Fcgr1", "Fcer1g", "Tyrobp", "Ctsb", "Lgmn", "Cd68"),
  "Immunoregulatory / negative feedback" = c("Il10", "Tgfb1", "Tnfaip3", "Dusp1", "Socs3", "Zfp36"),
  "Monocyte-to-macrophage transition" = c("Csf1r", "Mafb", "Adgre1", "C1qa", "C1qb", "C1qc", "Apoe"))

mac_genes <- list(
  "Inflammatory cytokines" = c("Il1b", "Tnf", "Il6", "Il12b", "Ptgs2", "Cd38"),
  "Inflammasome / pyroptosis" = c("Nlrp3", "Pycard", "Casp1", "Gsdmd"),
  "Interferon response" = c("Stat1", "Irf1", "Isg15", "Ifit1", "Ifit3", "Gbp2"),
  "Chemokine secretion / recruitment" = c("Ccl2", "Ccl3", "Ccl4", "Ccl7", "Cxcl2", "Cxcl9", "Cxcl10"),
  "Chemokine sensing / trafficking" = c("Ccr1", "Ccr2", "Ccr5", "Cx3cr1", "Cxcr4"),
  "Adhesion / motility" = c("Itgam", "Itgal", "Itgb2", "Cd44", "Lsp1", "Rac2"),
  "Antigen processing / presentation" = c("H2-Aa", "H2-Ab1", "H2-Eb1", "Cd74", "Ciita", "Ctss"),
  "Phagocytosis / lysosome" = c("Fcgr1", "Fcer1g", "Tyrobp", "Ctsb", "Lgmn", "Cd68"),
  "Efferocytosis / resolution" = c("Mertk", "Axl", "Timd4", "Lrp1", "Stab1", "Cd163"),
  "Immunoregulatory / negative feedback" = c("Il10", "Tgfb1", "Tnfaip3", "Dusp1", "Socs3", "Zfp36"),
  "Alternative activation / tissue repair" = c("Mrc1", "Arg1", "Retnla", "Chil3", "Egr2", "Igf1"),
  "ECM remodelling / angiogenesis" = c("Spp1", "Lgals3", "Mmp12", "Mmp14", "Fn1", "Vegfa"),
  "Lipid handling / homeostasis" = c("Trem2", "Apoe", "Lpl", "Abca1", "Fabp5", "Plin2"))

#categories in Figure 3; the rest go in Figure S6
mono_fig3 <- c("Inflammatory cytokines / alarmins", "Chemokine secretion / recruitment", "Chemokine sensing / trafficking",
               "Adhesion / transendothelial migration", "Patrolling / vascular surveillance",
               "Immunoregulatory / negative feedback", "Monocyte-to-macrophage transition")
mac_fig3 <- c("Inflammatory cytokines", "Chemokine secretion / recruitment", "Chemokine sensing / trafficking",
              "Antigen processing / presentation", "Phagocytosis / lysosome", "Efferocytosis / resolution",
              "Immunoregulatory / negative feedback", "Alternative activation / tissue repair",
              "Lipid handling / homeostasis")

sample_to_condition <- c(S1 = "Control", S5 = "Control", S3 = "PEPITEM", S4 = "PEPITEM", S7 = "PEPITEM", S8 = "PEPITEM")
samples <- names(sample_to_condition)

#mean log-normalised expression (counts per 10k, log1p) for a set of cells
mean_expr <- function(counts, cells) {
  m <- counts[, cells, drop = F]
  m <- m[, colSums(m) > 0, drop = F]
  m <- sweep(m, 2, colSums(m) / 10000, "/")
  list(avg = rowMeans(log1p(m)), pct = 100 * rowMeans(m > 0))
}

zscore <- function(x, lim) {
  f <- x[is.finite(x)]
  if (length(f) < 2 || sd(f) == 0) return(ifelse(is.na(x), NA, 0))
  pmax(-lim, pmin(lim, (x - mean(f)) / sd(f)))
}

```

```{r}

functional_panels <- function(obj, panel, fig3, types, name) {
  genes <- intersect(unlist(panel), rownames(obj))
  category <- setNames(rep(names(panel), lengths(panel)), unlist(panel))[genes]
  counts <- LayerData(obj[["RNA"]], layer = "counts")[genes, ]
  cond <- sample_to_condition[as.character(obj$Sample_all)]

  #dot colour: per condition within each subset (and all cells), z-scored per gene across everything, capped at 2
  groups <- expand.grid(Population = c("All", types), Condition = c("Control", "PEPITEM"), stringsAsFactors = F)
  dots <- bind_rows(lapply(seq_len(nrow(groups)), function(i) {
    cells <- colnames(obj)[cond == groups$Condition[i] & (groups$Population[i] == "All" | obj$population == groups$Population[i])]
    e <- mean_expr(counts, cells)
    data.frame(Gene = genes, Category = category, groups[i, ], Average = e$avg, Percent = e$pct)
  }))
  dots <- dots %>% group_by(Gene) %>% mutate(Scaled = zscore(Average, 2)) %>% ungroup()

  for (type in types) {
    #per-mouse mean, z-scored per gene across mice (capped at 1.5) -> heatmaps (Figs S4/S5)
    per_mouse <- sapply(samples, function(s) {
      cells <- colnames(obj)[obj$Sample_all == s & obj$population == type]
      if (length(cells)) mean_expr(counts, cells)$avg else rep(NA, length(genes))
    })
    z <- t(apply(per_mouse, 1, zscore, lim = 1.5))
    hm <- data.frame(Gene = rep(genes, length(samples)), Mouse = rep(samples, each = length(genes)), z = as.vector(z))
    hm$Gene <- factor(hm$Gene, levels = rev(genes))
    hm$Mouse <- factor(hm$Mouse, levels = samples)
    p <- ggplot(hm, aes(Mouse, Gene, fill = z)) + geom_tile() +
      scale_fill_gradient2(low = "#4575B4", mid = "#FFFFBF", high = "#D73027", limits = c(-1.5, 1.5)) +
      labs(title = type, x = NULL, y = NULL) + theme_minimal()
    ggsave(file.path(out, paste0("heatmap_", name, "_", gsub("[^A-Za-z]+", "_", type), ".pdf")), p, width = 4, height = 12)

    #Figure 3 genes: change = |mean z PEPITEM - mean z Control|; top 3 per category, only if change >= 0.8
    change <- abs(rowMeans(z[, sample_to_condition[colnames(z)] == "PEPITEM", drop = F], na.rm = T) -
                  rowMeans(z[, sample_to_condition[colnames(z)] == "Control", drop = F], na.rm = T))
    top <- unlist(lapply(fig3, function(cat) {
      g <- genes[category == cat]
      g <- g[order(-change[g])][seq_len(min(3, length(g)))]
      g[!is.na(change[g]) & change[g] >= 0.8]
    }))
    extra <- genes[!category %in% fig3]

    for (set in list(list(g = top, tag = "fig3"), list(g = extra, tag = "figS6"))) {
      d <- dots[dots$Population == type & dots$Gene %in% set$g, ]
      if (!nrow(d)) next
      d$Gene <- factor(d$Gene, levels = genes[genes %in% set$g])
      d$Category <- factor(d$Category, levels = names(panel))
      p <- ggplot(d, aes(Gene, Condition)) +
        geom_point(aes(size = Percent, fill = Scaled), shape = 21, colour = "grey25") +
        scale_fill_gradient2(low = "#0057D9", mid = "#FFFFFF", high = "#E31A1C", limits = c(-2, 2)) +
        scale_size(range = c(0.4, 4), limits = c(0, 100)) +
        facet_grid(. ~ Category, scales = "free_x", space = "free_x") +
        labs(title = type, x = NULL, y = NULL) + theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"), strip.text = element_text(size = 6))
      ggsave(file.path(out, paste0(set$tag, "_", name, "_", gsub("[^A-Za-z]+", "_", type), ".pdf")), p,
             width = 2 + 0.25 * length(set$g), height = 2.5)
    }
  }
  dots
}

mono_dots <- functional_panels(mono, mono_genes, mono_fig3, mono_types, "monocytes")
mac_dots <- functional_panels(mac, mac_genes, mac_fig3, mac_types, "macrophages")
write.csv(rbind(mono_dots, mac_dots), file.path(out, "functional_dotplot_values.csv"), row.names = F)

```

