#composition per mouse and propeller-style statistics
#needs cell_fate_per_cell.csv from subclustering.R

```{r}

library(dplyr)
library(limma)
library(ggplot2)

dir <- "/rds/projects/m/mcgetthm-rnaseq-arthritis/AU_final_analysis"
out <- file.path(dir, "AU_final_figures_update", "composition_stats")
dir.create(out, recursive = T, showWarnings = F)

fate <- read.csv(file.path(dir, "AU_final_figures_update", "exclusion_flow_sensitivity", "cell_fate_per_cell.csv"))
sample_to_condition <- c(S1 = "Control", S5 = "Control", S3 = "PEPITEM", S4 = "PEPITEM", S7 = "PEPITEM", S8 = "PEPITEM")
fate$Condition <- sample_to_condition[fate$Mouse]

mono_types <- c("Classical", "Transitional", "Non-classical")
mac_types <- c("MHC-II-high", "TLF+ resident-like")
fate$Broad <- fate$Final_population
fate$Broad[fate$Final_population %in% mono_types] <- "Monocytes"
fate$Broad[fate$Final_population %in% mac_types] <- "Macrophages"

#cells per mouse per population, retained cells only
count_cells <- function(d, col) d %>% count(Sample = Mouse, Condition, Cluster = .data[[col]], name = "Cluster_cells")
kept <- fate[fate$Status == "Retained", ]
total <- table(kept$Mouse)

broad <- count_cells(kept, "Broad")
subsets <- count_cells(kept[kept$Final_population %in% c(mono_types, mac_types), ], "Final_population")
subsets$Cluster <- recode(subsets$Cluster, "MHC-II-high" = "MHC-II-hi interstitial", "TLF+ resident-like" = "TLF+ resident")
broad$Total_cells <- as.numeric(total[broad$Sample])
subsets$Total_cells <- as.numeric(total[subsets$Sample])

sum(total)   # 22108

```

```{r}

#arcsine square-root proportions, limma with robust eBayes, BH within each family
test_props <- function(d) {
  d$Condition <- factor(d$Condition, levels = c("Control", "PEPITEM"))
  samp <- d %>% distinct(Sample, Condition, Total_cells) %>% arrange(Condition, Sample)
  clus <- unique(d$Cluster)
  m <- sapply(samp$Sample, function(s) sapply(clus, function(cl) {
    v <- d$Cluster_cells[d$Sample == s & d$Cluster == cl]; if (length(v)) v[1] else 0 }))
  m <- matrix(m, nrow = length(clus), dimnames = list(clus, samp$Sample))
  prop <- sweep(m, 2, samp$Total_cells, "/")
  fit <- eBayes(lmFit(asin(sqrt(prop)), model.matrix(~samp$Condition)), robust = T)
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  ctl <- samp$Condition == "Control"
  data.frame(Population = clus,
             Mean_percent_control = 100 * rowMeans(prop[, ctl, drop = F]),
             Mean_percent_PEPITEM = 100 * rowMeans(prop[, !ctl, drop = F]),
             P_value = tt$P.Value, Adjusted_P_value_BH = tt$adj.P.Val)
}

stats <- rbind(cbind(Family = "Broad populations", test_props(broad)),
               cbind(Family = "Monocyte and macrophage subsets", test_props(subsets)))
write.csv(stats, file.path(out, "composition_statistics_propeller.csv"), row.names = F)
stats

```

```{r}

#sensitivity - all recovered cells, with all excluded cells (removed neutrophils included) pooled into one "Excluded" group
fate$Broad_all <- ifelse(fate$Status == "Excluded", "Excluded", fate$Broad)
all_total <- table(fate$Mouse)

broad2 <- count_cells(fate, "Broad_all")
broad2$Total_cells <- as.numeric(all_total[broad2$Sample])
subsets2 <- subsets
subsets2$Total_cells <- as.numeric(all_total[subsets2$Sample])

sum(all_total)   # 32953

stats2 <- rbind(cbind(Family = "Broad populations", test_props(broad2)),
                cbind(Family = "Monocyte and macrophage subsets", test_props(subsets2)))
write.csv(stats2, file.path(out, "composition_statistics_propeller_all_cells.csv"), row.names = F)
stats2

```

```{r}

#bar = mean, error bar = mean + SD, points = individual mice
plot_props <- function(d) {
  d$Percent <- 100 * d$Cluster_cells / d$Total_cells
  s <- d %>% group_by(Cluster, Condition) %>% summarise(Mean = mean(Percent), SD = sd(Percent), .groups = "drop")
  ggplot(s, aes(Cluster, Mean, fill = Condition, colour = Condition)) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    geom_errorbar(aes(ymin = Mean, ymax = Mean + SD), position = position_dodge(0.8), width = 0.25) +
    geom_point(data = d, aes(Cluster, Percent), position = position_dodge(0.8), shape = 21, fill = NA) +
    scale_fill_manual(values = c(Control = "#AEE0FF", PEPITEM = "#FBC7BE")) +
    scale_colour_manual(values = c(Control = "#005C8C", PEPITEM = "#B3382C")) +
    labs(x = NULL, y = "Proportion of cells (%)") + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

ggsave(file.path(out, "composition_broad.pdf"), plot_props(broad), width = 6, height = 4)
ggsave(file.path(out, "composition_subsets.pdf"), plot_props(subsets), width = 5, height = 4)
ggsave(file.path(out, "composition_broad_all_cells.pdf"), plot_props(broad2), width = 6, height = 4)

```