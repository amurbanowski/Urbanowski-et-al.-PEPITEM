#pseudobulk DE (DESeq2) and GO for monocytes and macrophages
#needs cell_fate_per_cell.csv from subclustering.R

```{r}

library(Seurat)
library(Matrix)
library(DESeq2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(dplyr)
library(ggplot2)
library(ggrepel)

dir <- "/rds/projects/m/mcgetthm-rnaseq-arthritis/AU_final_analysis"
out <- file.path(dir, "AU_pseudo_revised")
dir.create(out, recursive = T, showWarnings = F)

fate <- read.csv(file.path(dir, "AU_final_figures_update", "exclusion_flow_sensitivity",
                           "cell_fate_per_cell.csv"))
fate <- fate[fate$Status == "Retained", ]

mono_types <- c("Classical", "Transitional", "Non-classical")
mac_types <- c("MHC-II-high", "TLF+ resident-like")

aggr <- readRDS(file.path(dir, "AU_input_objects", "all_synovial_object_used.rds"))
DefaultAssay(aggr) <- "RNA"
aggr <- aggr[, fate$Barcode[fate$Final_population %in% c(mono_types, mac_types)]]
if (inherits(aggr[["RNA"]], "Assay5")) aggr <- JoinLayers(aggr)
aggr$subtype <- fate$Final_population[match(colnames(aggr), fate$Barcode)]

sample_to_condition <- c(S1 = "Control", S5 = "Control", S3 = "PEPITEM", S4 = "PEPITEM", S7 = "PEPITEM", S8 = "PEPITEM")

table(aggr$subtype, aggr$Sample_all)

```

```{r}

#sum counts per mouse; mice with < 10 cells dropped; need 2+ mice per group
pseudo_de <- function(obj) {
  n <- table(as.character(obj$Sample_all))
  keep <- names(n)[n >= 10]
  cond <- factor(sample_to_condition[keep], levels = c("Control", "PEPITEM"))
  if (any(table(cond) < 2)) return(NULL)

  counts <- LayerData(obj[["RNA"]], layer = "counts")
  m <- sapply(keep, function(s) rowSums(counts[, obj$Sample_all == s, drop = F]))
  m <- round(m)
  m <- m[rowSums(m) >= 10, ]

  dds <- DESeqDataSetFromMatrix(m, data.frame(row.names = keep, Condition = cond), ~Condition)
  dds <- DESeq(dds, quiet = T)
  res <- results(dds, contrast = c("Condition", "PEPITEM", "Control"), alpha = 0.05)
  res <- lfcShrink(dds, contrast = c("Condition", "PEPITEM", "Control"), res = res, type = "normal")

  de <- data.frame(Gene = rownames(res), log2FC = res$log2FoldChange, P_value = res$pvalue,
                   Adjusted_P_value = res$padj)
  de$Direction <- "Not Sig"
  de$Direction[!is.na(de$Adjusted_P_value) & de$Adjusted_P_value <= 0.05 & de$log2FC >= 0.25] <- "Up"
  de$Direction[!is.na(de$Adjusted_P_value) & de$Adjusted_P_value <= 0.05 & de$log2FC <= -0.25] <- "Down"
  de[order(de$P_value), ]
}

comparisons <- list(
  All_monocytes = mono_types, Classical_monocytes = "Classical", Transitional_monocytes = "Transitional",
  Non_classical_monocytes = "Non-classical",
  All_macrophages = mac_types, MHC_II_high_macrophages = "MHC-II-high",
  TLF_resident_like_macrophages = "TLF+ resident-like")

de <- list()
for (cmp in names(comparisons)) {
  de[[cmp]] <- pseudo_de(aggr[, aggr$subtype %in% comparisons[[cmp]]])
  write.csv(de[[cmp]], file.path(out, paste0("DE_", cmp, ".csv")), row.names = F)
}

#DEGs (Down + Up): All monocytes 5, Classical 2, Transitional 0, Non-classical 0,
#All macrophages 10, MHC-II-high 7, TLF+ resident-like 5
sapply(de, function(x) table(factor(x$Direction, levels = c("Down", "Up"))))

```

```{r}

#GO - genes at nominal P < 0.05 (no fold-change cut), up and down separately
#background = all genes tested; BP, CC and MF; terms of 5-500 genes; BH
to_entrez <- function(x) {
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (!length(x)) return(character())
  m <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db, keys = x, columns = "ENTREZID", keytype = "SYMBOL"))
  unique(na.omit(m$ENTREZID))
}

odds_ratio <- function(k, n, M, N) {
  v <- c(k, n - k, M - k, N - M - (n - k))
  if (any(v < 0)) return(NA)
  if (any(v == 0)) v <- v + 0.5
  v[1] * v[4] / (v[2] * v[3])
}

run_go <- function(d) {
  universe <- to_entrez(d$Gene[is.finite(d$P_value)])
  bind_rows(lapply(c(DOWN = -1, UP = 1), function(sgn) {
    genes <- intersect(to_entrez(d$Gene[is.finite(d$P_value) & d$P_value < 0.05 & sign(d$log2FC) == sgn]), universe)
    if (length(genes) < 3) return(NULL)
    t <- as.data.frame(enrichGO(genes, universe = universe, OrgDb = org.Mm.eg.db, keyType = "ENTREZID",
                                ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                                minGSSize = 5, maxGSSize = 500))
    t <- t[!is.na(t$p.adjust) & t$p.adjust <= 0.05, ]
    if (!nrow(t)) return(NULL)
    r <- function(x) as.numeric(sub("/.*", "", x)); tot <- function(x) as.numeric(sub(".*/", "", x))
    t$Odds_ratio <- mapply(odds_ratio, r(t$GeneRatio), tot(t$GeneRatio), r(t$BgRatio), tot(t$BgRatio))
    t$Direction <- if (sgn < 0) "DOWN" else "UP"
    t[is.finite(t$Odds_ratio), ]
  }))
}

go <- list()
for (cmp in names(de)) {
  if (is.null(de[[cmp]])) next
  go[[cmp]] <- run_go(de[[cmp]])
  if (!is.null(go[[cmp]]) && nrow(go[[cmp]])) write.csv(go[[cmp]], file.path(out, paste0("GO_", cmp, ".csv")), row.names = F)
}

sapply(go, function(x) if (is.null(x)) 0 else nrow(x))

```

```{r}

#terms shown in the figures: BP only, inflammation/migration terms first then by adjusted P,
#overlapping terms (Jaccard > 0.5 on annotated genes) collapsed into the higher one, 10 per direction
#overlap uses all tested genes, as in the published figures (enrichGO itself uses genes with a P value)
focus <- "inflamm|cytokine|interferon|chemokine|myeloid activation|monocyte activation|macrophage activation|leukocyte activation|chemotaxis|migration|motility|adhesion|extravasation|transendothelial|phagocyt|efferocyt|apoptotic cell clearance|wound healing|resolution|anti-inflammatory"

go_sets <- function(ids, universe) {
  t <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db, keys = unique(ids), keytype = "GOALL", columns = "SYMBOL"))
  t <- t[!is.na(t$SYMBOL) & t$SYMBOL %in% universe, ]
  split(t$SYMBOL, t$GOALL)
}

pick_terms <- function(g, universe, n = 10) {
  g <- g[g$ONTOLOGY == "BP", ]
  if (!nrow(g)) return(g)
  g <- g[order(g$Direction, !grepl(focus, g$Description, ignore.case = T), g$p.adjust, -g$Odds_ratio), ]
  sets <- go_sets(g$ID, universe)
  bind_rows(lapply(split(g, g$Direction), function(d) {
    kept <- character()
    for (id in d$ID) {
      gi <- sets[[id]]
      if (is.null(gi) || !any(vapply(kept, function(k) length(intersect(gi, sets[[k]])) / length(union(gi, sets[[k]])) > 0.5, logical(1))))
        kept <- c(kept, id)
    }
    head(d[d$ID %in% kept, ], n)
  }))
}

go_plot <- function(g, title) {
  g$Description <- factor(g$Description, levels = rev(unique(g$Description)))
  ggplot(g, aes(Direction, Description)) +
    geom_point(aes(size = Count, fill = Odds_ratio), shape = 21, colour = "grey25") +
    scale_fill_gradientn(colours = c("#FFF200", "#FFB000", "#FF6A00", "#EF2B1F", "#B2182B"), trans = "log10") +
    labs(x = NULL, y = NULL, title = title, fill = "Odds ratio", size = "Genes") + theme_bw()
}

volcano <- function(d, title) {
  d <- d[!is.na(d$Adjusted_P_value), ]
  lab <- d[d$Direction != "Not Sig", ]
  ggplot(d, aes(log2FC, -log10(Adjusted_P_value), colour = Direction)) +
    geom_point(size = 0.8) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed") +
    geom_text_repel(data = lab, aes(label = Gene), colour = "black", fontface = "italic", max.overlaps = Inf) +
    scale_colour_manual(values = c(Down = "#2C7FB8", "Not Sig" = "#BDBDBD", Up = "#B2182B")) +
    labs(x = "Average log2 fold change", y = "-log10 adjusted P", title = title) + theme_bw()
}

for (cmp in names(de)) {
  if (is.null(de[[cmp]])) next
  ggsave(file.path(out, paste0("volcano_", cmp, ".pdf")), volcano(de[[cmp]], cmp), width = 5, height = 5)
  if (!is.null(go[[cmp]]) && nrow(go[[cmp]])) {
    shown <- pick_terms(go[[cmp]], de[[cmp]]$Gene)
    if (nrow(shown)) ggsave(file.path(out, paste0("GO_", cmp, ".pdf")), go_plot(shown, cmp), width = 7, height = 5)
  }
}

```
