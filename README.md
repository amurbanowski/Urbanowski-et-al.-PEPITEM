# scRNA-seq

Code for the higher-resolution re-analysis of synovial monocytes and macrophages from PEPITEM-treated mice with antigen-induced arthritis.

Data are from Wahid et al. (GEO: GSE324735). Sample processing, QC, doublet removal, integration (RPCA) and broad clustering/annotation were done with `process_samples.R` and `analysis.R` from https://github.com/chrismahony/Wahid-et-al.-PEPITEM. The integrated object from `analysis.R` (`all_synovial_object_used.rds`) is the starting point here.

Run in this order:

1. `process_samples.R` - QC, doublet removal, integration (RPCA) and clustering (from: https://github.com/chrismahony/Wahid-et-al.-PEPITEM)
2. `analysis.R` - broad cell annotation (from: https://github.com/chrismahony/Wahid-et-al.-PEPITEM)
3. `subclustering.R` - monocyte and macrophage subclustering, then the final label for every cell
4. `composition.R` - composition per mouse and stats
5. `pseudobulk.R` - pseudobulk DE, GO, volcano and GO plots
6. `figures.R` - UMAPs, marker dot plots and functional gene panels

Paths are set for our HPC, change `dir` at the top of each script.
