library('DESeq2')
library('dplyr')
library('stringr')
library('tibble')
library("EnhancedVolcano")
library("pheatmap")
library(pasilla)
library("vsn")
library(fuzzyjoin)
library("factoextra")
library("NbClust")
library(data.table)
library("GenomicFeatures")
library("purrr")
library("cluster")
library("gprofiler2")
library('clusterProfiler')
library('enrichplot')
library('org.Hs.eg.db')

LOADING_TO_DESEQ2 <- function(merged_qorts_input, pre_filtering = FALSE) {
  
  # Read data from indir and start reformatting to DESeq2 format
  dds_DF = read.delim(merged_qorts_input, header = FALSE, sep = "\t")

  if (ncol(dds_DF) == 10) {
    DF_colnames = c("gene_id", "NES1", "NES2", "NES3", "Progenitor1", "Progenitor2", "Progenitor3", "Neuron1", "Neuron2", "Neuron3")
  } else if (ncol(dds_DF) == 7) {
    DF_colnames = c("gene_id", "NES1", "NES2", "NES3", "Neuron1", "Neuron2", "Neuron3")
  }
  
  colnames(dds_DF) = DF_colnames
  dds_DF = tibble::column_to_rownames(dds_DF, "gene_id")
  
  # Load annotation file and save gene id and gene name for downstream conversion of gene id to gene name
  # Get column data. lapply creates lists of the column names, i.e list of cell types and replicate number
  # The list are used do produce the data frame coldata where replicates are nested under the cell types.
  # The rownames of coldata are the column names of All_RNA_DF
  
  if (ncol(dds_DF) == 9) {
    groups            = lapply(colnames(dds_DF), function(x) {substr(x, 1, nchar(x)-1)}) %>% unlist %>% factor(levels = c("NES", "Progenitor", "Neuron"))
    replicates        = lapply(colnames(dds_DF), function(x) {substr(x, nchar(x), nchar(x)+1)}) %>% unlist
    coldata           = data.frame(groups = as.factor(groups), replicates = as.factor(replicates))
    rownames(coldata) = colnames(dds_DF)
  } else if (ncol(dds_DF) == 6) {
    groups            = lapply(colnames(dds_DF), function(x) {substr(x, 1, nchar(x)-1)}) %>% unlist %>% factor(levels = c("NES", "Neuron"))
    replicates        = lapply(colnames(dds_DF), function(x) {substr(x, nchar(x), nchar(x)+1)}) %>% unlist
    coldata           = data.frame(groups = as.factor(groups), replicates = as.factor(replicates))
    rownames(coldata) = colnames(dds_DF)
  }
  
  # Load data into DESeq experiement object
  # Load the DESeq data set by providing the countData which is All_RNA_DF. The output data is organised by the input matrix which is coldata.
  # Design expresses how the counts for each gene depend on the variable in coldata.
  dds = DESeqDataSetFromMatrix(countData = dds_DF,
                               colData = coldata,
                               design = ~ groups)
  
  # Depending on RNA library size, prefiltering of low read counts can be applied. See below for more information.
  # FROM DESeq2 vignette:
  # "While it is not necessary to pre-filter low count genes before running the DESeq2 functions, there are two reasons which make pre-filtering useful: 
  # by removing rows in which there are very few reads, we reduce the memory size of the dds data object, and we increase the speed of the transformation 
  # and testing functions within DESeq2."
  # Prefiltering step:
  
  if (pre_filtering == TRUE) {
    keep = rowSums(counts(dds))>=10
    dds = dds[keep,]
  }
  
  # Calculate DESeq transformation and obtain results, FDR correction
  dds = DESeq(dds)
  
  return(dds)
  
}

CONVERT_ENSEMBL_TO_GENE_ID <- function(input_matrix = assay(sdf), annotation_file, 
                                       subset_genes = selected_genes, create_DF = FALSE) {
  
  annotation_load_file = readRDS(annotation_file)
  annotation_load_file = as.data.frame(annotation_load_file) %>% dplyr::select(., gene_id, gene_name)
  annotation_load_file$gene_id = gsub("\\..*","", annotation_load_file$gene_id)
  
  if (is.matrix(input_matrix) == TRUE) {
    output_DF = as.data.frame(input_matrix) %>% tibble::rownames_to_column(., "row") %>% 
      tibble::add_column(., gene_name = 0, .before = 1)
    output_DF$row = gsub("\\..*","", output_DF$row)
    output_DF$gene_name = annotation_load_file$gene_name[match(output_DF$row, annotation_load_file$gene_id)]
    output_DF = subset(output_DF, select = -c(row))
  } else if (has_rownames(input_matrix) == FALSE & is.matrix(input_matrix) == FALSE) {
    output_DF = as.data.frame(input_matrix) %>% rename(., row = gene_name) %>% tibble::add_column(., gene_name = 0, .before = 1)
    output_DF$gene_name = annotation_load_file$gene_name[match(output_DF$row, annotation_load_file$gene_id)]
    output_DF = subset(output_DF, select = -c(row))
  }
  
  if (is.null(subset_genes) == FALSE) {
    output_DF = filter(output_DF, gene_name %in% subset_genes)
  }
  output_DF = output_DF %>% tibble::remove_rownames(.) %>% tibble::column_to_rownames(., var="gene_name")

  if (create_DF == TRUE){
    return(output_DF)
  } else if (create_DF == FALSE) {
    output_matrix = as.matrix(output_DF)
    return(output_matrix)
  }
}

NORMALISATION_TPM = function(dds_count = counts(dds)) {
  
  TxDb = makeTxDbFromUCSC(genome = "hg19", tablename = "ensGene")
  TxDb = sum(width(GenomicRanges::reduce(exonsBy(TxDb, by="gene"))))
  TxDb = as.data.table(TxDb, keep.rownames = TRUE)
  TxDb$TxDb = TxDb$TxDb/1000
  
  count_TPM = as.data.table(dds_count, keep.rownames = TRUE)
  count_TPM$rn = gsub("\\..*","", count_TPM$rn)
  count_TPM = left_join(count_TPM, TxDb, by="rn")
  
  ### CHECK IF THE MAPPLY PART IS WORKING ###

  count_TPM[,2:(ncol(count_TPM)-1)] = count_TPM[,2:(ncol(count_TPM)-1)] / count_TPM[,ncol(count_TPM)]
  scale_factor = (base::colSums(count_TPM[2:(ncol(count_TPM)-1)], na.rm = TRUE))/1000000
  count_TPM[,2:(ncol(count_TPM)-1)] = mapply('/', count_TPM[,2:(ncol(count_TPM)-1)], scale_factor)
  count_TPM[,ncol(count_TPM)] = NULL
  rownames(count_TPM) = count_TPM[,1]
  count_TPM[,1] <- NULL

  return(count_TPM)
  
}

FILTER_TURNED_GENES <- function(res_dds_list = res, alpha = 0.05, log2FC_limit = 2, strict_filter = FALSE, 
                                tpm_count = count_norm, tpm_limit = 1,  delete_item = NULL) {
  
  if (is.null(delete_item) == FALSE){
    res[-c(delete_item)]
  }
  
  select_list = list()
  print(paste0("Filtering with Log2 FC limit: ", log2FC_limit))
  
  if (length(res_dds_list) == 3) {
    if (strict_filter == FALSE & is.null(tpm_count) == FALSE) {
      
      tpm_count = rownames_to_column(tpm_count, var = "rn")

      select_list[["NES_up"]] = which(res[["NsPr"]]$padj <= alpha & res[["NsPr"]]$log2FoldChange < -log2FC_limit & 
                                        res[["PrNe"]]$padj <= alpha & 
                                        res[["NsNe"]]$log2FoldChange < -log2FC_limit &
                                        tpm_count[,1] >= tpm_limit & tpm_count[,2] >= tpm_limit & tpm_count[,3] >= tpm_limit)
      select_list[["Progenitor_up"]] = which(res[["NsPr"]]$padj <= alpha & res[["NsPr"]]$log2FoldChange > log2FC_limit & 
                                               res[["PrNe"]]$padj <= alpha & res[["PrNe"]]$log2FoldChange < -log2FC_limit & 
                                               tpm_count[,4] >= tpm_limit & tpm_count[,5] >= tpm_limit & tpm_count[,6] >= tpm_limit)
      select_list[["Neuron_up"]] = which(res[["PrNe"]]$padj <= alpha & res[["PrNe"]]$log2FoldChange > log2FC_limit & 
                                           res[["NsPr"]]$padj <= alpha & 
                                           res[["NsNe"]]$log2FoldChange > log2FC_limit &
                                           tpm_count[,7] >= tpm_limit & tpm_count[,8] >= tpm_limit & tpm_count[,9] >= tpm_limit)
      select_list[["NES_down"]] = which(res[["NsPr"]]$padj <= alpha & res[["NsPr"]]$log2FoldChange > log2FC_limit & 
                                          res[["PrNe"]]$padj <= alpha & 
                                          res[["NsNe"]]$log2FoldChange > log2FC_limit &
                                          tpm_count[,1] >= tpm_limit & tpm_count[,2] >= tpm_limit & tpm_count[,3] >= tpm_limit)
      select_list[["Progenitor_down"]] = which(res[["NsPr"]]$padj <= alpha & res[["NsPr"]]$log2FoldChange < -log2FC_limit & 
                                                 res[["PrNe"]]$padj <= alpha & res[["PrNe"]]$log2FoldChange > log2FC_limit & 
                                                 tpm_count[,4] >= tpm_limit & tpm_count[,5] >= tpm_limit & tpm_count[,6] >= tpm_limit)
      select_list[["Neuron_down"]] = which(res[["PrNe"]]$padj <= alpha & res[["PrNe"]]$log2FoldChange < -log2FC_limit & 
                                             res[["NsPr"]]$padj <= alpha & 
                                             res[["NsNe"]]$log2FoldChange < -log2FC_limit &
                                             tpm_count[,7] >= tpm_limit & tpm_count[,8] >= tpm_limit & tpm_count[,9] >= tpm_limit)
      
    } else if (strict_filter == TRUE) {
      tpm_count = rownames_to_column(tpm_count, var = "rn")
      
      select_list[["NES_up"]] = which(res[["NsPr"]]$padj <= alpha & res[["NsPr"]]$log2FoldChange < -log2FC_limit & 
                                        res[["PrNe"]]$padj <= alpha & res[["PrNe"]]$log2FoldChange > -log2FC_limit & res[["PrNe"]]$log2FoldChange < log2FC_limit &
                                        res[["NsNe"]]$log2FoldChange < -log2FC_limit & res[["NsNe"]]$padj <= alpha &
                                        tpm_count[,1] >= tpm_limit & tpm_count[,2] >= tpm_limit & tpm_count[,3] >= tpm_limit)
      select_list[["Progenitor_up"]] = which(res[["NsPr"]]$padj <= alpha & res[["NsPr"]]$log2FoldChange > log2FC_limit & 
                                               res[["PrNe"]]$padj <= alpha & res[["PrNe"]]$log2FoldChange < -log2FC_limit &
                                               res[["NsNe"]]$padj <= alpha & res[["NsNe"]]$log2FoldChange > -log2FC_limit & res[["NsNe"]]$log2FoldChange < log2FC_limit &
                                               tpm_count[,4] >= tpm_limit & tpm_count[,5] >= tpm_limit & tpm_count[,6] >= tpm_limit)
      select_list[["Neuron_up"]] = which(res[["PrNe"]]$padj <= alpha & res[["PrNe"]]$log2FoldChange > log2FC_limit & 
                                           res[["NsPr"]]$padj <= alpha & res[["NsPr"]]$log2FoldChange > -log2FC_limit & res[["NsPr"]]$log2FoldChange < log2FC_limit & 
                                           res[["NsNe"]]$log2FoldChange > log2FC_limit & res[["NsNe"]]$padj <= alpha &
                                           tpm_count[,7] >= tpm_limit & tpm_count[,8] >= tpm_limit & tpm_count[,9] >= tpm_limit)
      select_list[["NES_down"]] = which(res[["NsPr"]]$padj <= alpha & res[["NsPr"]]$log2FoldChange > log2FC_limit & 
                                          res[["PrNe"]]$padj <= alpha & res[["PrNe"]]$log2FoldChange > -log2FC_limit & res[["PrNe"]]$log2FoldChange < log2FC_limit &
                                          res[["NsNe"]]$log2FoldChange > log2FC_limit & res[["NsNe"]]$padj <= alpha &
                                          tpm_count[,1] >= tpm_limit & tpm_count[,2] >= tpm_limit & tpm_count[,3] >= tpm_limit)
      select_list[["Progenitor_down"]] = which(res[["NsPr"]]$padj <= alpha & res[["NsPr"]]$log2FoldChange < -log2FC_limit & 
                                                 res[["PrNe"]]$padj <= alpha & res[["PrNe"]]$log2FoldChange > log2FC_limit &
                                                 res[["NsNe"]]$padj <= alpha & res[["NsNe"]]$log2FoldChange > -log2FC_limit & res[["NsNe"]]$log2FoldChange < log2FC_limit &
                                                 tpm_count[,4] >= tpm_limit & tpm_count[,5] >= tpm_limit & tpm_count[,6] >= tpm_limit)
      select_list[["Neuron_down"]] = which(res[["PrNe"]]$padj <= alpha & res[["PrNe"]]$log2FoldChange < -log2FC_limit & 
                                             res[["NsPr"]]$padj <= alpha & res[["NsPr"]]$log2FoldChange > -log2FC_limit & res[["NsPr"]]$log2FoldChange < log2FC_limit & 
                                             res[["NsNe"]]$log2FoldChange < -log2FC_limit & res[["NsNe"]]$padj <= alpha &
                                             tpm_count[,7] >= tpm_limit & tpm_count[,8] >= tpm_limit & tpm_count[,9] >= tpm_limit)
      
    }
  } else if (length(res_dds_list) == 1) {
    if (strict_filter == FALSE & is.null(tpm_count) == FALSE) {
      
      res[["NsNe"]] = as.data.frame(res$NsNe) %>% rownames_to_column(., var = "rn")
      tpm_count = rownames_to_column(tpm_count, var = "rn")
      res[["NsNe"]] = res$NsNe %>% right_join(tpm_count, ., by = "rn") %>% column_to_rownames(., var = "rn")
      
      select_list[["NES_up"]] = which(res[["NsNe"]]$padj <= alpha & res[["NsNe"]]$log2FoldChange < -log2FC_limit & res[["NsNe"]][,1] >= tpm_limit &
                                        res[["NsNe"]][,2] >= tpm_limit & res[["NsNe"]][,3] >= tpm_limit)
      select_list[["Neuron_up"]] = which(res[["NsNe"]]$padj <= alpha & res[["NsNe"]]$log2FoldChange > log2FC_limit & res[["NsNe"]][,4] >= tpm_limit &
                                           res[["NsNe"]][,5] >= tpm_limit & res[["NsNe"]][,6] >= tpm_limit)
      select_list[["no_change"]] = which(res[["NsNe"]]$padj <= alpha & res[["NsNe"]]$log2FoldChange < log2FC_limit & res[["NsNe"]]$log2FoldChange > -log2FC_limit & 
                                           res[["NsNe"]][,1] >= tpm_limit & res[["NsNe"]][,2] >= tpm_limit & res[["NsNe"]][,3] >= tpm_limit &
                                           res[["NsNe"]][,4] >= tpm_limit & res[["NsNe"]][,5] >= tpm_limit & res[["NsNe"]][,6] >= tpm_limit)
      
      
      
    } else if (strict_filter == TRUE) {
      select_list[["NES_up"]] = which(res[["NsPr"]]$padj <= alpha & res[["NsPr"]]$lfcSE < -log2FC_limit & res[["PrNe"]]$padj >= alpha)
      select_list[["Neuron_up"]] = which(res[["PrNe"]]$padj <= alpha & res[["PrNe"]]$lfcSE > log2FC_limit & res[["NsPr"]]$padj >= alpha)
      select_list[["NES_down"]] = which(res[["NsPr"]]$padj <= alpha & res[["NsPr"]]$lfcSE > log2FC_limit & res[["PrNe"]]$padj >= alpha)
      select_list[["Neuron_down"]] = which(res[["PrNe"]]$padj <= alpha & res[["PrNe"]]$lfcSE < -log2FC_limit & res[["NsPr"]]$padj > alpha)
    }
  }
  
  
  return(select_list)
  
}
TARGET_GENE_LIST <- function(gene_targets, annotation_file = annotation_file, 
                             sdf = sdf) {

  output_DF = as.data.frame(assay(sdf)) %>% tibble::rownames_to_column(., "row") %>% 
    tibble::add_column(., gene_name = 0, .before = 1) 
  output_DF$gene_name = annotation_file$gene_name[match(output_DF$row, annotation_file$gene_id)]
  output_DF = subset(output_DF, select = -c(row))
  output_DF = output_DF %>% filter(., gene_name %in% gene_targets)
  output_DF = output_DF %>% tibble::remove_rownames(.) %>% tibble::column_to_rownames(., var="gene_name")
  output_matrix = as.matrix(output_DF)
  
  return(gene_DF)
  
}

PLOT_FILTERED_TURNED_GENES <- function(select_list = NULL,
                                     filter_regulation = FALSE,
                                     selected_data_transformation = sdf,
                                     pheatmap_annotation = dds_selected_df,
                                     annotation = annotation_file, 
                                     show_gene_names = FALSE, 
                                     cluster_cols = FALSE, 
                                     log2FC_suffix) {
  plot_list = list()
  for (reg in 1:length(select_list)) {
    
      if (length(select_list[[reg]]) < 2) {
        next
      }
    
      select_sdf = assay(sdf)[select_list[[reg]],]
      reg_suffix = str_extract(names(select_list[reg]), "^[:alpha:]*")
      type_suffix = str_extract(names(select_list[reg]), "[:alpha:]*$")
      
      if (show_gene_names == FALSE) {
      plot_list[[names(select_list[reg])]] = pheatmap(select_sdf, cluster_rows = TRUE, show_rownames = FALSE, 
                                                      cluster_cols = cluster_cols, annotation_col = pheatmap_annotation, 
                                                      main = paste0(type_suffix, " expressed genes in ", reg_suffix,
                                                                    "\n Log2 fold change threshold = ", log2FC_suffix))
      } else if (show_gene_names == TRUE) {
        
        select_converted_sdf = CONVERT_ENSEMBL_TO_GENE_ID(input_matrix = select_sdf, annotation_file = annotation, 
                                                   create_DF = FALSE, subset_genes = NULL)
        
        plot_list[[names(select_list[reg])]] = pheatmap(select_converted_sdf, cluster_rows = TRUE, show_rownames = TRUE, 
                                                        cluster_cols = cluster_cols, annotation_col = pheatmap_annotation, 
                                                        main = paste0(type_suffix, " expressed genes in ", reg_suffix,
                                                                      "\n Log2 fold change threshold = ", log2FC_suffix))
      }
        
    } 
  return(plot_list)
}

PLOT_SELECTED_GENES <- function(selected_genes_sdf, pheatmap_annotation, cluster_cols = FALSE) {
  
  pheatmap_plot = pheatmap(selected_genes_sdf, cluster_rows = TRUE, show_rownames = TRUE, 
                           cluster_cols = cluster_cols, annotation_col = pheatmap_annotation, 
                           main = paste0("Gene regulation of selected genes across cell types"))
  
  return(pheatmap_plot)
  
}

SAMPLE_DIST_HEATMAP <- function(df) {
  
  sampleDists <- dist(t(assay(df)))
  sampleDistMatrix <- as.matrix( sampleDists )
  library( "gplots" )
  library( "RColorBrewer" )
  colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  mat_cluster_rows = hclust(sampleDists)
  mat_cluster_rows = dendextend::rotate(mat_cluster_rows, c(2,3))

  mat_cluster_cols = hclust(sampleDists)
  mat_cluster_cols = dendextend::rotate(mat_cluster_cols, c(2,3))



  heatmap_plot = pheatmap( sampleDistMatrix, trace="none", col=colours, 
                           main = "Heatmap of sample distances", 
                           cluster_rows = mat_cluster_rows,
                           cluster_cols = mat_cluster_cols)

  return(heatmap_plot)
  
}

K_CLUSTER_OPTIMISATION = function(select_filtered_list, data_trans = sdf, indir, annotation_file) {
  
  set.seed(123)
  
  for (i in 1:length(select_filtered_list)) {
    
    df = assay(data_trans)[select_filtered_list[[i]], ]
    df = CONVERT_ENSEMBL_TO_GENE_ID(input_matrix = df, create_DF = FALSE, subset_genes = NULL, annotation_file = annotation_file)
    
    distance = get_dist(df, method = "euclidean")
    distance_plot = fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
    ggsave(filename = paste0(names(select_filtered_list[i]), "_euclidean_distance_heatmap.png"), plot = distance_plot, device = "png", path = "QC/")
    
    if (nrow(df) <= 2) {
      next
    }
    
    # Elbow method
    elbow_plot = fviz_nbclust(df, kmeans, method = "wss")
    
    # Silhouette method
    silhouette_plot = fviz_nbclust(df, kmeans, method = "silhouette")
    
    # Gap statistic method
    gap_stat <- clusGap(df, FUN = kmeans, nstart = 25,
                        K.max = 10, B = 50)
    gap_kmean_plot = fviz_gap_stat(gap_stat)
    
    gap_stat <- clusGap(df, FUN = hcut, nstart = 25,
                        K.max = 10, B = 50)
    gap_hc_plot = fviz_gap_stat(gap_stat)
    
    ggsave(filename = paste0(names(select_filtered_list[i]), "_k_clusters_elbow.png"), plot = elbow_plot, device = "png", path = "QC/")
    ggsave(filename = paste0(names(select_filtered_list[i]), "_k_clusters_silhouette.png"), plot = silhouette_plot, device = "png", path = "QC/")
    ggsave(filename = paste0(names(select_filtered_list[i]), "_k_clusters_gap_kmean.png"), plot = gap_kmean_plot, device = "png", path = "QC/")
    ggsave(filename = paste0(names(select_filtered_list[i]), "_k_clusters_gap_hc.png"), plot = gap_hc_plot, device = "png", path = "QC/")
    
  }
  
}

RNA_clusterProfiler_selected = function(res_dds, select_filter, tpm_filter, filter_suffix, output_dir, 
                                        alpha = 0.05, logthreshold = 0.5, basemean_threshold = 10) {
  
  universe_df = as.data.frame(res_dds, keep.rownames = TRUE) %>% rownames_to_column(., var = "gene_id")
  universe_df$gene_id = gsub("\\..*","", universe_df$gene_id)
  
  ## feature 1: numeric vector, i.e. log2 FC
  geneList = universe_df[,3]
  
  ## feature 2: named vector, i.e. ensembl name
  names(geneList)  = as.character(universe_df[,1])
  
  ## Subsetting filtered genes
  gene = geneList[select_filter]
  
  ## Filtering universe for TPM >= 1
  geneList = geneList[tpm_filter]
  
  ## feature 3: Sorting both lists to decreasing orders
  geneList = sort(geneList, decreasing = TRUE)
  gene = sort(gene, decreasing = TRUE)
  
  ## Running GO clssification
  geneList.df = clusterProfiler::bitr(names(geneList), fromType = "ENSEMBL", 
                                      toType = c("ENTREZID", "SYMBOL"),
                                      OrgDb = "org.Hs.eg.db")
  
  gene.df = clusterProfiler::bitr(names(gene), fromType = "ENSEMBL", 
                                  toType = c("ENTREZID", "SYMBOL"),
                                  OrgDb = "org.Hs.eg.db")
  
  ggo_CC <- groupGO(gene     = names(gene),
                    OrgDb    = org.Hs.eg.db,
                    keyType = "ENSEMBL",
                    ont      = "CC",
                    level    = 7,
                    readable = TRUE)
  
  ggo_MF <- groupGO(gene     = names(gene),
                    OrgDb    = org.Hs.eg.db,
                    keyType = "ENSEMBL",
                    ont      = "MF",
                    level    = 3,
                    readable = TRUE)
  
  ggo_BP <- groupGO(gene     = names(gene),
                    OrgDb    = org.Hs.eg.db,
                    keyType = "ENSEMBL",
                    ont      = "BP",
                    level    = 4,
                    readable = TRUE)
  
  ggo_CC_plot = barplot(ggo_CC, drop=TRUE, showCategory=10, order = TRUE) + 
    labs(title = paste0("GO classification of cellular compartments in ", filter_suffix, "\n n = ", length(gene)),
         y ="Gene Count")
  ggo_MF_plot = barplot(ggo_MF, drop=TRUE, showCategory=10, order = TRUE) + 
    labs(title = paste0("GO classification of molecular function in ", filter_suffix, "\n n = ", length(gene)),
         y ="Gene Count")
  ggo_BP_plot = barplot(ggo_BP, drop=TRUE, showCategory=10, order = TRUE) + 
    labs(title = paste0("GO classification of biological processes in ", filter_suffix, "\n n = ", length(gene)),
         y ="Gene Count")
  
  ggsave(filename = paste0(filter_suffix, "_ggo_CC_plot.pdf"), 
         plot = ggo_CC_plot, device = "pdf", path = output_dir,
         limitsize = FALSE, dpi = "retina", width = 10, height = 8)
  ggsave(filename = paste0(filter_suffix, "_ggo_MF_plot.pdf"), 
         plot = ggo_MF_plot, device = "pdf", path = output_dir,
         limitsize = FALSE, dpi = "retina", width = 10, height = 8)
  ggsave(filename = paste0(filter_suffix, "_ggo_BP_plot.pdf"), 
         plot = ggo_BP_plot, device = "pdf", path = output_dir,
         limitsize = FALSE, dpi = "retina", width = 10, height = 8)
  
  # Running enrich GO
  ego_CC <- enrichGO(gene       = names(gene),
                     universe      = names(geneList),
                     OrgDb         = org.Hs.eg.db,
                     keyType       = "ENSEMBL",
                     ont           = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
  
  ego_MF <- enrichGO(gene          = names(gene),
                     universe      = names(geneList),
                     OrgDb         = org.Hs.eg.db,
                     keyType       = "ENSEMBL",
                     ont           = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
  
  ego_BP <- enrichGO(gene          = names(gene),
                     universe      = names(geneList),
                     OrgDb         = org.Hs.eg.db,
                     keyType       = "ENSEMBL",
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
  
  ego_CC_plot = barplot(ego_CC, showCategory=10, order = TRUE) + 
    labs(title = paste0("Enriched GO of cellular compartments in ", filter_suffix, "\n n = ", length(gene)),
         y ="Gene Count") 
  
  ego_MF_plot = barplot(ego_MF, showCategory=10) + 
    labs(title = paste0("Enriched GO of molecular functions in ", filter_suffix, "\n n = ", length(gene)),
         y ="Gene Count")
  
  ego_BP_plot = barplot(ego_BP, showCategory=10) + 
    labs(title = paste0("Enriched GO of biological processes in ", filter_suffix, "\n n = ", length(gene)),
         y ="Gene Count")
  
  ggsave(filename = paste0(filter_suffix, "_ego_CC_plot.pdf"), 
         plot = ego_CC_plot, device = "pdf", path = output_dir,
         limitsize = FALSE, dpi = "retina", width = 10, height = 8)
  ggsave(filename = paste0(filter_suffix, "_ego_MF_plot.pdf"), 
         plot = ego_MF_plot, device = "pdf", path = output_dir,
         limitsize = FALSE, dpi = "retina", width = 10, height = 8)
  ggsave(filename = paste0(filter_suffix, "_ego_BP_plot.pdf"), 
         plot = ego_BP_plot, device = "pdf", path = output_dir,
         limitsize = FALSE, dpi = "retina", width = 10, height = 8)
  
  ego_list = list()
  if ((nrow(ego_CC) > 0) == TRUE) {
    ego_CC_df = as.data.frame(ego_CC)
    ego_list[["CC"]] = ego_CC_df
  }
  if ((nrow(ego_MF) > 0) == TRUE) {
    ego_MF_df = as.data.frame(ego_MF)
    ego_list[["MF"]] = ego_MF_df
  }
  if ((nrow(ego_BP) > 0) == TRUE) {
    ego_BP_df = as.data.frame(ego_BP)
    ego_list[["BP"]] = ego_BP_df
  }
  
  return(ego_list)
  
}

# WORKS
RNA_clusterProfiler_gseGO = function(res_dds, tpm_filter, suffix_comp, output_dir) {
  
  universe_df = as.data.frame(res_dds, keep.rownames = TRUE) %>% rownames_to_column(., var = "gene_id") #%>%
  universe_df$gene_id = gsub("\\..*","", universe_df$gene_id)
  
  ## feature 1: numeric vector, i.e. log2 FC
  geneList = universe_df[,3]
  
  ## feature 2: named vector, i.e. ensembl name
  names(geneList)  = as.character(universe_df[,1])
  
  ## Filtering universe for TPM >= 1
  geneList = geneList[tpm_filter]
  
  ## feature 3: Sorting both lists to decreasing orders
  geneList = sort(geneList, decreasing = TRUE)
  
  gse_CC <- gseGO(geneList     = geneList,
                  OrgDb        = org.Hs.eg.db,
                  keyType       = "ENSEMBL",
                  ont          = "CC",
                  nPerm        = 1000,
                  minGSSize    = 100,
                  maxGSSize    = 500,
                  pvalueCutoff = 0.05,
                  verbose      = FALSE)
  
  gse_MF <- gseGO(geneList     = geneList,
                  OrgDb        = org.Hs.eg.db,
                  keyType       = "ENSEMBL",
                  ont          = "MF",
                  nPerm        = 1000,
                  minGSSize    = 100,
                  maxGSSize    = 500,
                  pvalueCutoff = 0.05,
                  verbose      = FALSE)
  
  gse_BP <- gseGO(geneList     = geneList,
                  OrgDb        = org.Hs.eg.db,
                  keyType       = "ENSEMBL",
                  ont          = "BP",
                  nPerm        = 1000,
                  minGSSize    = 100,
                  maxGSSize    = 500,
                  pvalueCutoff = 0.05,
                  verbose      = FALSE)
  
  gse_CC_plot = clusterProfiler::dotplot(gse_CC, showCategory = 10, title = paste0("GSEA of GO regarding cellular compartments in ", suffix_comp))
  gse_MF_plot = clusterProfiler::dotplot(gse_MF, showCategory = 10, title = paste0("GSEA of GO regarding molecular functions in ", suffix_comp))
  gse_BP_plot = clusterProfiler::dotplot(gse_BP, showCategory = 10, title = paste0("GSEA of GO regarding biological processes in ", suffix_comp))
  
  ggsave(filename = paste0(suffix_comp, "_gse_CC_plot.pdf"), 
         plot = gse_CC_plot, device = "pdf", path = output_dir,
         limitsize = FALSE, dpi = "retina", width = 10, height = 8)
  ggsave(filename = paste0(suffix_comp, "_gse_MF_plot.pdf"), 
         plot = gse_MF_plot, device = "pdf", path = output_dir,
         limitsize = FALSE, dpi = "retina", width = 10, height = 8)
  ggsave(filename = paste0(suffix_comp, "_gse_BP_plot.pdf"), 
         plot = gse_BP_plot, device = "pdf", path = output_dir,
         limitsize = FALSE, dpi = "retina", width = 10, height = 8)
  
  gse_list = list()
  if (nrow(gse_CC) > 0) {
    gse_CC_df = as.data.frame(gse_CC)
    gse_list[["CC"]] = gse_CC_df
  }
  if (nrow(gse_MF) > 0) {
    gse_MF_df = as.data.frame(gse_MF)
    gse_list[["MF"]] = gse_MF_df
  }
  if (nrow(gse_BP) > 0) {
    gse_BP_df = as.data.frame(gse_BP)
    gse_list[["BP"]] = gse_BP_df
  }
  
  return(gse_list)
  
}

