# Description: RNA analysis script which filters RNA data for different expression
# patterns across cell types.

indir = "/RNA-BLISS_Data_analysis_Pipeline/Output/"
annotation_file = "/RNA-BLISS_Data_analysis_Pipeline/Data/hg19/annotation/hg19_Gencode19_annotations.all.genes.regions.rds"
merged_qorts_dir = "/RNA-BLISS_Data_analysis_Pipeline/Data/RNA_data/RNA_qorts/qorts_formatted/B193/All_merged_RNA.qorts.txt"
save_plots = TRUE
log2FC_threshold = 1
tpm = 1
alpha = 0.05
Apply_strict = TRUE
plot_QC = TRUE
save_plot_gene_id = FALSE
k_mean_analysis = FALSE

  
library('cowplot')
library('ggplot2')
library('ggpubr')
library("pheatmap")
library("dendextend")
library('IHW')

source("/RNA-BLISS_Data_analysis_Pipeline/Scripts/RNA_analysis/RNA_analysis_tools.R")

dir.create(paste0(indir, "RNA_analysis/"), showWarnings = FALSE)
dir.create(paste0(indir, "RNA_analysis/LFC_", log2FC_threshold, "_TPM_", tpm), showWarnings = FALSE)
indir = paste0(indir, "RNA_analysis/LFC_", log2FC_threshold, "_TPM_", tpm)
setwd(indir)
dir.create("results", showWarnings = FALSE)
dir.create("results/k-mean_clustering", showWarnings = FALSE)
dir.create("results/DE_Heatmaps", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)
dir.create("results/plots", showWarnings = FALSE)
dir.create("results/GSEA", showWarnings = FALSE)
#dir.create(c("results/k-mean_clustering", "results/DE_Heatmaps", "results/tables", "results/plots"))

#strict_filter = TRUE
dds = LOADING_TO_DESEQ2(merged_qorts_input = merged_qorts_dir, 
                        pre_filtering = FALSE)

# Generating normalised count table
count_norm = NORMALISATION_TPM(dds_count = counts(dds))
fwrite(count_norm, file = paste0("tables/RNA_TPM_normalised_counts.csv"), sep = "\t", row.names = TRUE)
saveRDS(count_norm, file = paste0("tables/RNA_TPM_normalised_counts.rds"))

if (ncol(count_norm) == 6) {
  
  tpm_filter = which(count_norm[,1] >= 1 & count_norm[,2] >= 1 & count_norm[,3] >= 1 & 
                       count_norm[,4] >= 1 & count_norm[,5] >= 1 & count_norm[,6] >= 1)
  
} else if (ncol(count_norm) == 9) {
  
  tpm_filter = which(count_norm[,1] >= 1 & count_norm[,2] >= 1 & count_norm[,3] >= 1 & 
                       count_norm[,4] >= 1 & count_norm[,5] >= 1 & count_norm[,6] >= 1 &
                       count_norm[,7] >= 1 & count_norm[,8] >= 1 & count_norm[,9] >= 1)
  
}

#Setting up for loop
res = list()

# Independent hypothesis weighting to optimize power
resIHW = results(dds, filterFun = ihw)
summary(resIHW)
sum(resIHW$pvalue < alpha, na.rm = TRUE)

# Apending each cell type with replicates to res list. Log2 fold change and shrinked log2 fold change estimated 
# is determined. The ashr distribution is applied
res[["NsNe"]] = results(dds, alpha = alpha, contrast=c("groups", "Neuron", "NES"), pAdjustMethod = "fdr")
res[["NsNe"]] = lfcShrink(dds, contrast=c("groups", "Neuron", "NES"), res=res[["NsNe"]], type = "ashr")

if (nlevels(dds$groups) == 3) {
  res[["NsPr"]] = results(dds, alpha = alpha, contrast=c("groups", "Progenitor", "NES"), pAdjustMethod = "fdr")
  res[["NsPr"]] = lfcShrink(dds, contrast=c("groups", "Progenitor", "NES"), res=res[["NsPr"]], type = "ashr")
  res[["PrNe"]] = results(dds, alpha = alpha, contrast=c("groups", "Neuron", "Progenitor"), pAdjustMethod = "fdr")
  res[["PrNe"]] = lfcShrink(dds, contrast=c("groups", "Neuron", "Progenitor"), res=res[["PrNe"]], type = "ashr")
}

# Quality control, checking count data transformation
if (plot_QC == TRUE) {
  
  dir.create("QC", showWarnings = FALSE)
  
  # Producing boxplot of Cooks distance. Cook’s distance is a measure of how much a single sample is influencing 
  # the fitted coefficients for a gene, and a large value of Cook’s distance is intended to indicate an outlier count.
  #pdf(filename = "QC/Cooks_boxplot.pdf")
  #par(mar=c(8,5,2,2))
  #boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2, ylab="Cook's distance", main = "Cook's distance of samples")
  #dev.off()
  
  #melt_Cooks = melt(as.data.table(log10(assays(dds)[["cooks"]])))
  #ggplot(melt_Cooks, aes(x = variable, y = value)) + stat_boxplot(geom = "errorbar", width = 0.5) + geom_boxplot() + 
  #  labs(title = "Cook's distance of samples", y = "Cook's distance", x = "") +
    #main("Cook's distance of samples") + ylab("Cook's distance") + xlab("") +
  #  theme_cowplot() +
  #  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  #ggsave(filename = "Cooks_boxplot.pdf", plot = Cooks_boxplot, device = "pdf", path = "QC/")
  
  
  # SD should increase as the mean increase
  vsd <- vst(dds, blind=FALSE)
  rld <- rlog(dds, blind=FALSE)
  ntd <- normTransform(dds)
  ntd_plot = meanSdPlot(assay(ntd))
  rld_plot = meanSdPlot(assay(rld))
  vsd_plot = meanSdPlot(assay(vsd))
  ggsave(filename = "Sd_plot_ntd.pdf", plot = ntd_plot$gg, device = "pdf", path = "QC/")
  ggsave(filename = "Sd_plot_rld.pdf", plot = rld_plot$gg, device = "pdf", path = "QC/")
  ggsave(filename = "Sd_plot_vsd.pdf", plot = vsd_plot$gg, device = "pdf", path = "QC/")
  
  sampleDist_heat = SAMPLE_DIST_HEATMAP(df = vsd)
  ggsave(filename = "vsd_Euclidean_sample_dist.pdf", plot = sampleDist_heat, device = "pdf", path = "QC/")
  PCA_plot = plotPCA(vsd, intgroup = c("groups", "replicates")) + theme_cowplot()
  ggsave(filename = "vsd_PCA_plot.pdf", plot = PCA_plot, path = "QC/")
  
  sampleDist_heat = SAMPLE_DIST_HEATMAP(df = rld)
  ggsave(filename = "rld_Euclidean_sample_dist.pdf", plot = sampleDist_heat, device = "pdf", path = "QC/")
  PCA_plot = plotPCA(rld, intgroup = c("groups", "replicates")) + theme_cowplot()
  ggsave(filename = "rld_PCA_plot.pdf", plot = PCA_plot, path = "QC/")
  
  sampleDist_heat = SAMPLE_DIST_HEATMAP(df = ntd)
  ggsave(filename = "ntd_Euclidean_sample_dist.pdf", plot = sampleDist_heat, device = "pdf", path = "QC/")
  PCA_plot = plotPCA(ntd, intgroup = c("groups", "replicates")) + theme_cowplot()
  ggsave(filename = "ntd_PCA_plot.pdf", plot = PCA_plot, path = "QC/")
  
  sdf = vsd #sdf = selected_data_transformation
  
  
  drawLines <- function() abline(h=c(-log2FC_threshold,log2FC_threshold),col="dodgerblue",lwd=2)
  
  if (length(res) > 1) {
    for (i in 1:length(res)) {
      png(filename = paste0("QC/", names(res[i]), "_MA_plot.png"))
      plotMA(res[[i]], ylim=c(-10,10), alpha = 0.05, main = paste(names(res[i])), ylab = "Log2 fold change", xlab = "Mean of normalised counts");drawLines()
      dev.off()
    }
  } 
  
  #png(filename = "QC/MA_plots.png")
  pdf(file = "QC/MA_plots.pdf")
  if (length(res) > 1) {
    par(mfrow=c(1,3), oma = c(0, 0, 2, 0))
    for (i in 1:length(res)) {
      plotMA(res[[i]], ylim=c(-10,10), alpha = 0.05, main = paste(names(res[i])), ylab = "Log2 fold change", xlab = "Mean of normalised counts");drawLines()

    }
  } else if (length(res) == 1) {
    plotMA(res[[1]], ylim=c(-10,10), main = paste(names(res[1])));drawLines()
  }
  title(main = paste0("MA plot of log fold change across normalised means \n Log2 FC threshold = ", log2FC_threshold), outer = TRUE)
  dev.off()
  
  ggsave("MA_plot_single.pdf", plot = MA_plot, device = "pdf", path = "QC/")
  
  for (i in 1:length(res)) {
    
    png(filename = paste0("QC/",names(res[i]),"_p-value_normalised-mean.png"))
    plot(res[[i]]$baseMean+1, -log10(res[[i]]$pvalue),
         log="x", xlab="mean of normalized counts",
         ylab=expression(-log[10](pvalue)),
         ylim=c(0,30),
         cex=.4, col=rgb(0,0,0,.3),
         main = paste(names(res[i])))
    dev.off()
  }
  
} else if (plot_QC == FALSE) {
  vsd <- vst(dds, blind=FALSE)
  # Default is vsd
  sdf = vsd
}


## Check which has the best data transformation.
# For current RNA data, vsd is used.
dds_selected_df = as.data.frame(colData(dds)[,c("groups", "replicates")])
#selected_genes = c("RBFOX3", "NES", "SOX2", "TBR1", "EOMES", "PAX6", "MAP2", "DLG4", "SYN1", "GAD65", "pcr1", "PCR2", 
#                   "EZH1", "EZH2", "BRCA1", "BRCA2", "MYCN", "MYC", "TOP2A", "TOP2B", "UBE3A", "MLL", "GAPDH", "KMT2A",
#                   "DCX", "XIST", "WWOX", "EMX1", "MECP2", "JARID1C", "HUWE1", "RLIM", "POL2")

selected_genes = c("DAB1", "NEGR1", "LPHN2", "PRKG1", "PCDH15", "CTNNA3", "NRG3", "NAV2", "LRRC4C", "DLG2", "SOX5", "GPC6",
                   "MDGA2", "RBFOX1", "DCC", "CTNNA2", "NCKAP5", "LRP1B", "ERBB4", "MACROD2", "LARGE", "ERC2", "LSAMP", "LPP",
                   "CTNND2", "CDH18", "PARK2", "SDK1", "AUTS2", "MAGI2", "EXOC4", "CSMD1", "LINGO2", "MID1", "IL1RAPL1", "PCDH11X")
select_filtered_list = FILTER_TURNED_GENES(res_dds_list = res, alpha = alpha, log2FC_limit = log2FC_threshold, strict_filter = Apply_strict,
                                           tpm_count = count_norm, tpm_limit = tpm, delete_item = NULL)

plot_select_filtered_list = PLOT_FILTERED_TURNED_GENES(select_list = select_filtered_list,
                                                       filter_regulation = TRUE,
                                                       selected_data_transformation = sdf,
                                                       annotation = annotation_file,
                                                       pheatmap_annotation = dds_selected_df,
                                                       show_gene_names = FALSE,
                                                       cluster_cols = FALSE,
                                                       log2FC_suffix = log2FC_threshold)

plot_select_filtered_names_list = PLOT_FILTERED_TURNED_GENES(select_list = select_filtered_list,
                                                             filter_regulation = TRUE,
                                                             selected_data_transformation = sdf,
                                                             annotation = annotation_file,
                                                             pheatmap_annotation = dds_selected_df,
                                                             show_gene_names = TRUE,
                                                             cluster_cols = FALSE,
                                                             log2FC_suffix = log2FC_threshold)

# Gene Set Enrichment Analysis section
GSEA_dir = paste0(getwd(), "/results/GSEA/")
# Running GSEA GO using all genes in comparisions
for (i in 1:length(res)) {

  gseGO_list = RNA_clusterProfiler_gseGO(res_dds = res[[i]], 
                                         tpm_filter = tpm_filter,
                                         suffix_comp = names(res[i]),
                                         output_dir = GSEA_dir)
  
  for (j in 1:length(gseGO_list)) {
    saveRDS(gseGO_list[[j]], file = paste0(GSEA_dir, names(res[i]), "_gseGO_",  names(gseGO_list[j]), "_table.rds"))
  }
  
}

# Running GSEA enrich GO and group GO, requiring the filtered out genes.
# For the analysis, a minimum of 10 selected genes are required. 
for (i in 1:length(select_filtered_list)) {
  
  if ((length(select_filtered_list[[i]]) > 10) == TRUE) {
    enrichGO_list = RNA_clusterProfiler_selected(res_dds = res$NsNe, 
                                                 select_filter = select_filtered_list[[i]], 
                                                 tpm_filter = tpm_filter,
                                                 filter_suffix = names(select_filtered_list[i]), 
                                                 output_dir = GSEA_dir)
    
    for (j in 1:length(enrichGO_list)) {
      saveRDS(enrichGO_list[[j]], file = paste0(GSEA_dir, "enrichGO_",  names(enrichGO_list[j]), "_table.rds"))
    }
  } 
  
}

if (k_mean_analysis == TRUE) {
  
  # K-mean clustering section
  k_opt = K_CLUSTER_OPTIMISATION(select_filtered_list = select_filtered_list, data_trans = sdf, indir = indir, annotation_file = annotation_file)
  
  # k_clusters has to have as many values as length of selected_filtered_list
  k_clusters_v = c(3, 0, 3, 3, 3, 3)
  k_means_res_list = list()
  set.seed(123)
  for (i in 1:length(select_filtered_list)) {
    
    if (length(select_filtered_list[[i]]) <= 2) {
      next
    }
    
    count_filt = assay(sdf)[select_filtered_list[[i]],]
    
    count_filt_label = CONVERT_ENSEMBL_TO_GENE_ID(input_matrix = count_filt, create_DF = FALSE, subset_genes = NULL, annotation_file = annotation_file)
    
    k_mean_res = kmeans(count_filt, k_clusters_v[[i]], nstart = 25)
    k_mean_res_clusters = as.data.table(k_mean_res$cluster, keep.rownames = "rn") %>% dplyr::rename(., rn = V1) %>% dplyr::rename(., k_cluster = V2)
    
    k_mean_res = kmeans(count_filt_label, k_clusters_v[[i]], nstart = 25)
    
    k_mean_plot = fviz_cluster(k_mean_res, data = count_filt, geom = "point", main = paste0(names(select_filtered_list[i]), " k-means clustering plot")) +
      theme_cowplot()
    ggsave(filename = paste0(names(select_filtered_list[i]), "_k-means_clustering_plot.pdf"), plot = k_mean_plot, device = "pdf", path = "results/k-mean_clustering/")
    
    k_mean_label_plot = fviz_cluster(k_mean_res, data = count_filt, main = paste0(names(select_filtered_list[i]), " k-means clustering plot"))  +
      theme_cowplot()
    ggsave(filename = paste0(names(select_filtered_list[i]), "_k-means_clustering_labels_plo.pdf"), plot = k_mean_label_plot, device = "pdf", path = "results/k-mean_clustering/")
    
  }
  
}

selected_genes = CONVERT_ENSEMBL_TO_GENE_ID(input_matrix = assay(sdf), annotation_file = annotation_file, 
                                            subset_genes = selected_genes, create_DF = FALSE)

plot_selected_genes = PLOT_SELECTED_GENES(selected_genes_sdf = selected_genes, 
                                          pheatmap_annotation = dds_selected_df,
                                          cluster_cols = FALSE)


if (is.null(plot_selected_genes) == FALSE) {
  ggsave(filename = "selected_genes_pheatmap.pdf", plot = plot_selected_genes, device = "pdf", path = "results/DE_Heatmaps/")
}

save_plot_gene_id = FALSE
#i = 1
barplot_list = list()
gene_count_lineplot_list = list()

if (save_plots == TRUE) {
  
  for (i in 1:length(select_filtered_list)) {
    
    if (length(select_filtered_list[[i]]) < 2) {
      next
    }
    
    suffix_i = (names(select_filtered_list[i]))
    
    counts_i = as.data.table(counts(dds, normalized = TRUE)[select_filtered_list[[i]], ], keep.rownames = "gene_name")
    
    if (save_plot_gene_id == TRUE) {
      
      if (ncol(counts_i) == 10) {
        counts_i = as.data.table(CONVERT_ENSEMBL_TO_GENE_ID(input_matrix = counts_i, annotation_file = annotation_file), keep.rownames = "gene_name")
        counts_i = counts_i %>% mutate(mean_NES = round(rowMeans(.[,2:4]), digits = 2)) %>% mutate(mean_Progenitor = round(rowMeans(.[,5:7]), digits = 2)) %>% 
          mutate(mean_Neuron = round(rowMeans(.[,8:10]), digits = 2)) %>% mutate(basemean = round(rowMeans(.[,2:10]), digits = 2))
      } else if (ncols(counts_i) == 7) {
        counts_i = as.data.table(CONVERT_ENSEMBL_TO_GENE_ID(input_matrix = counts_i, annotation_file = annotation_file), keep.rownames = "gene_name")
        counts_i = counts_i %>% mutate(mean_NES = round(rowMeans(.[,2:4]), digits = 2)) %>% mutate(mean_Neuron = round(rowMeans(.[,5:7]), digits = 2)) %>% 
          mutate(basemean = round(rowMeans(.[,2:7]), digits = 2))
      }
      
      fwrite(counts_i, file = paste0("tables/", suffix_i, "_dds_counts_gene_id.csv"), sep = "\t")

    } else if (save_plot_gene_id == FALSE) {
      
      counts_i$gene_name = gsub("\\..*","", counts_i$gene_name)
      
      if (ncol(counts_i) == 10) {
        counts_i = counts_i %>% mutate(mean_NES = round(rowMeans(.[,2:4]), digits = 2)) %>% mutate(mean_Progenitor = round(rowMeans(.[,5:7]), digits = 2)) %>% 
          mutate(mean_Neuron = round(rowMeans(.[,8:10]), digits = 2)) %>% mutate(basemean = round(rowMeans(.[,2:10]), digits = 2))
      } else if (ncols(counts_i) == 7) {
        counts_i = counts_i %>% mutate(mean_NES = round(rowMeans(.[,2:4]), digits = 2)) %>% mutate(mean_Neuron = round(rowMeans(.[,5:7]), digits = 2)) %>% 
          mutate(basemean = round(rowMeans(.[,2:7]), digits = 2))
      }
      fwrite(counts_i, file = paste0("tables/",suffix_i, "_dds_raw_counts.csv"), sep = "\t")
      saveRDS(counts_i, file = paste0("tables/",suffix_i, "_dds_raw_counts.rds"))
    }
    
    if (ncol(counts_i) == 14) {
      line_plot_genes_dt = counts_i[-c(2:10, 14)]
      colnames(line_plot_genes_dt) = c("Gene", "NES", "Progenitor", "Neuron")
      n_rows_line_plot_genes_dt = nrow(line_plot_genes_dt)
      line_plot_genes_dt = melt(line_plot_genes_dt)
      
      line_plot_gene = ggplot(line_plot_genes_dt, aes(x = variable, y = value, col = Gene, group = Gene)) +
        geom_point() + geom_line(size = 1) + theme_cowplot() + theme(legend.position = "none") +
        labs(title = paste0(suffix_i, " mean normalised read count per gene \n across all cell types. n = ", n_rows_line_plot_genes_dt), x="Cell type", 
             y = "Mean normalised read count")
      ggsave(filename = paste0(suffix_i, "_count_gene_lineplot.pdf"), plot = line_plot_gene, 
             device = "pdf", path = "results/plots/")
      gene_count_lineplot_list[[i]] = line_plot_gene
      
      line_plot_counts = c(sum(counts_i$mean_NES), sum(counts_i$mean_Progenitor), sum(counts_i$mean_Neuron))
      line_plot_names = c("NES", "Progenitor", "Neuron")
      line_plot_df = as.data.frame(cbind(line_plot_names, line_plot_counts))
      line_plot_df$line_plot_names = factor(line_plot_df$line_plot_names, 
                                            levels = line_plot_df$line_plot_names)
      line_plot_df$line_plot_counts = as.numeric(as.character(line_plot_df$line_plot_counts))

      read_count_plot = ggplot(line_plot_df, aes(x = as.factor(line_plot_names), y=line_plot_counts)) +
        geom_bar(stat = "identity", colour = "black") + theme_cowplot() + 
        labs(title = paste0("Sum of ",suffix_i, " mean read count"), x="Cell type", y = "Sum of mean read count")
      ggsave(filename = paste0(suffix_i, "_read_count_barplot.pdf"), plot = read_count_plot, 
             device = "pdf", path = "results/plots/")
      barplot_list[[i]] = read_count_plot
    } else if (ncol(counts_i) == 10) {
      line_plot_genes_dt = as.data.frame(counts_i[-c(2:7, 10)])
      colnames(line_plot_genes_dt) = c("Gene", "NES", "Neuron")
      line_plot_genes_dt = melt(line_plot_genes_dt)
      
      line_plot_gene = ggplot(line_plot_genes_dt, aes(x = variable, y = value, col = Gene, group = Gene)) +
        geom_point() + geom_line(size = 1) + theme_cowplot() + theme(legend.position = "none") +
        labs(title = paste0(suffix_i, " mean read count per gene \n across all cell types"), x="Cell type", 
             y = "Gene mean read count")
      ggsave(filename = paste0(suffix_i, "_count_gene_barplot.pdf"), plot = line_plot_gene, 
             device = "pdf", path = "results/plots/")
      gene_count_lineplot_list[[i]] = line_plot_gene
      
      line_plot_counts = c(sum(counts_i$mean_NES), sum(counts_i$mean_Neuron))
      line_plot_names = c("NES", "Neuron")
      line_plot_df = as.data.frame(cbind(line_plot_names, line_plot_counts))
      line_plot_df$line_plot_names = factor(line_plot_df$line_plot_names, 
                                            levels = line_plot_df$line_plot_names)
      line_plot_df$line_plot_counts = as.numeric(as.character(line_plot_df$line_plot_counts))

      read_count_plot = ggplot(line_plot_df, aes(x = as.factor(line_plot_names), y=line_plot_counts)) +
        geom_bar(stat = "identity", colour = "black") + theme_cowplot() + 
        labs(title = paste0("Sum of ",suffix_i, " mean read count"), x="Cell type", y = "Sum of mean read count")
      ggsave(filename = paste0(suffix_i, "_read_count_barplot.pdf"), plot = read_count_plot, 
             device = "pdf", path = "results/plots/")
      barplot_list[[i]] = read_count_plot
    }
    
  }
  
  if (ncol(counts_i) == 14) {
    arrange_bar_plot_count = ggarrange(plotlist = barplot_list, ncol = 3, nrow = 2)
    ggsave(filename = "All_cells_count_bar_plot.pdf", plot = arrange_bar_plot_count, 
           device = "pdf", width = 20, height = 15)
    arrange_line_plot_gene_count = ggarrange(plotlist = gene_count_lineplot_list, ncol = 3, nrow = 2)
    ggsave(filename = "All_cells_gene_count_line_plot.pdf", plot = arrange_line_plot_gene_count, 
           device = "pdf", width = 20, height = 15, path = "results/plots/")
  } else if (ncol(counts_i) == 10) {
    arrange_bar_plot_count = ggarrange(plotlist = barplot_list, ncol = 1, nrow = 2)
    ggsave(filename = "All_cells_count_bar_plot.pdf", plot = arrange_bar_plot_count, 
           device = "pdf", width = 20, height = 15, path = "results/plots/")
    arrange_line_plot_gene_count = ggarrange(plotlist = gene_count_lineplot_list, ncol = 1, nrow = 2)
    ggsave(filename = "All_cells_gene_count_line_plot.pdf", plot = arrange_line_plot_gene_count, 
           device = "pdf", width = 20, height = 15, path = "results/plots/")
  }
  
  for (i in 1:length(plot_select_filtered_list)) {
    
    ggsave(filename = paste0(names(plot_select_filtered_list[i]),"_pheatmap.pdf"), plot = plot_select_filtered_list[[i]], 
           device = "pdf", path = "results/DE_Heatmaps/")
    
  }
  
  for (i in 1:length(plot_select_filtered_names_list)) {
    
    ggsave(filename = paste0(names(plot_select_filtered_names_list[i]),"_pheatmap_names.pdf"), plot = plot_select_filtered_names_list[[i]], 
           device = "pdf", path = "results/DE_Heatmaps/")
    
  }
  
}
