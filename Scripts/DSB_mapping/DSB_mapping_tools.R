library("Rmisc")
library('data.table')
library('ggplot2')
library('ggpubr')
library('dplyr')
library('tidyr')
library('cowplot')
library('scales')
library("ggsignif")
library("ggpubr")
library("GenomicRanges")
library("rtracklayer")
library("VennDiagram")
library('Rsubread')

EDIT_COUNT_DT <- function(count_dt_dir) {
  
  dt_count = fread(file = count_dt_dir)
  dt_count = separate(dt_count, window, c("chr", "start", "end"), sep = "[ ]|[-]") %>%
    mutate(., Sum_Count = rowSums(.[,4:ncol(.)], na.rm = TRUE))
  dt_count$start = as.numeric(dt_count$start)
  dt_count$end = as.numeric(dt_count$end)
  dt_count$chr = gsub('chr', '', dt_count$chr)
  
  ### CURRENT SCRIPT REMOVES LAST BIN OF THE CHROMOSOME
  bin_size = dt_count$end[1] - dt_count$start[1]
  dt_count = subset(dt_count, (dt_count$end - dt_count$start) == bin_size)
  
  return(dt_count)
  
}

PLOT_DSB_CPM_ALL_chr <- function(CPM_DT, plot_title, window_size, BLISS_run) {
  
  CPM_DT$chr = factor(CPM_DT$chr, levels=unique(CPM_DT$chr))
  
  ggplot(CPM_DT) +
    geom_boxplot(aes(y = Sum_Count, x = chr)) +
    geom_hline(aes(yintercept = median(Sum_Count)), color = "red") +
    scale_x_discrete(limits = c(levels(CPM_DT$chr))) +
    ggtitle(paste0("Boxplots of CPM of DSB across ", plot_title, " chromosomes \nBin size = ", window_size)) +
    xlab("chromosome") + ylab("DSB count per million") + theme(plot.title = element_text(size = 15)) + theme_cowplot()
  
  ggsave(paste0(plot_title, "_CPM_All_chr_barplot.pdf"), width = 10, height = 7,  dpi = 500, path = paste0(BLISS_run, "/", window_size, "/",plot_title, "_chr_plots"))
  
  ggplot(CPM_DT) +
    geom_violin(aes(y = Sum_Count, x = chr)) +
    geom_hline(aes(yintercept = median(Sum_Count)), color = "red") +
    stat_summary(fun.y = "mean", geom="point", aes(y = Sum_Count, x = chr), color = "red") +
    stat_summary(fun.y = "median", geom="point", aes(y = Sum_Count, x = chr), color = "blue") +
    scale_x_discrete(limits = c(levels(CPM_DT$chr))) +
    ggtitle(paste0("Violin plot of CPM of DSB across ", plot_title, " chromosomes \nBin size = ", window_size)) +
    xlab("chromosome") + ylab("DSB count per million") + theme(plot.title = element_text(size = 15)) + theme_cowplot()
  
  ggsave(paste0(plot_title, "_CPM_All_chr_violinplot.pdf"), width = 10, height = 7,  dpi = 500, path = paste0(BLISS_run, "/", window_size, "/",plot_title, "_chr_plots"))
  
  
}

CALLING_BIN_GENE <- function(selected_bins, log2FC_threshold = log2(1.5), suffix_comp) {
  
  chr_vector = c(1:22, "X")
  
  selected_bins = subset(selected_bins, Log2fc >= log2FC_threshold | Log2fc <= -log2FC_threshold)
  selected_bins = as.data.table(selected_bins)
  #anno_chr = selected_bins$chr[1]
  
  anno = "/RNA-BLISS_Data_analysis_Pipeline/Data/hg19/gencode.v19.annotation.gtf"
  
  anno_DT = as.data.table(rtracklayer::import(anno))
  anno_DT$seqnames = gsub('chr', '', anno_DT$seqnames)
  #anno_DT = anno_DT[anno_DT$seqnames == anno_chr]
  anno_DT = anno_DT[anno_DT$type == "CDS"]
  anno_DT$gene_id = gsub("\\..*","", anno_DT$gene_id)
  anno_DT$transcript_id = gsub("\\..*","", anno_DT$transcript_id)
  #anno_DT = anno_DT[anno_DT$seqnames == chr_vector]
  
  setkey(selected_bins, chr, start, end)
  overlap_DT = foverlaps(anno_DT, selected_bins, by.x = c("seqnames", "start", "end"), type="any", which=TRUE)
  
  
  overlap_DT = foverlaps(selected_bins, anno_DT, by.x = start, by.y = end)
  #overlap_DT = overlap_DT[overlap_DT$yid != NA]
  overlap_DT = na.omit(overlap_DT)
  select_anno = anno_DT[overlap_DT$xid,]
  select_anno = select_anno[!duplicated(select_anno$gene_id)]
  #anno_DT = anno_DT[anno_DT$seqnames == gsub('chr', '', anno_DT$seqnames)]
  
}

RETURN_FEATURE_COUNTS_LIST <- function(bam_files, bam_indir, 
                        n_bam = NULL, keep_End = FALSE) {
  
  anno = "/RNA-BLISS_Data_analysis_Pipeline/Data/hg19/gencode.v19.annotation.gtf"
  
  fC_res_list = list()
  if (is.null(n_bam) == TRUE) {
    n_bam = length(bam_files)
  } 
  
  ## As long as we only have sBLISS for NES, the script is only used for NES and therefore only
  ## the first three BAM files are used
  #for (i in 1:length(bam_files)) {
  for (i in 1:n_bam) {
    
    fC_res  = featureCounts(files = paste0(bam_indir, bam_files[i]), annot.ext = anno, isGTFAnnotationFile = TRUE,
                            GTF.featureType = "exon", GTF.attrType = "transcript_id", GTF.attrType.extra = "gene_id")
    
    
    fC_res_dt = as.data.table(fC_res$annotation) %>% left_join(., as.data.table(fC_res$counts, keep.rownames = "GeneID"), by = "GeneID")
    names(fC_res_dt)[ncol(fC_res_dt)] = "Count"
    fC_res_dt = fC_res_dt[(fC_res_dt$Count > 0),]
    fC_res_dt = fC_res_dt[(grepl("chrM", fC_res_dt$Chr) == FALSE & grepl("chrY", fC_res_dt$Chr) == FALSE),]
    fC_res_dt$Length = fC_res_dt$Length/1000
    fC_res_dt$TPM = (fC_res_dt$Count/fC_res_dt$Length)
    sum_fC = sum(fC_res_dt$TPM)/1000000
    fC_res_dt$TPM = fC_res_dt$TPM/sum_fC
    
    if (keep_End == FALSE) {
      fC_res_dt = fC_res_dt[c(1:3, 7, 9)]
    } else if (keep_End == TRUE) {
      fC_res_dt = fC_res_dt[c(1:4, 7, 9)]
    }
    
    fC_res_dt = fC_res_dt[order(fC_res_dt$TPM, decreasing = TRUE),]
    fC_res_dt = fC_res_dt[(!duplicated(fC_res_dt$gene_id) == TRUE),]
    fC_res_dt = fC_res_dt[fC_res_dt$TPM >= 1,]
    
    if (n_bam == 6 | n_bam == 3) {
      
      if (i <= 3) {
        fC_res_list[[paste0("NES_", i)]] = fC_res_dt
      } else if (i > 3) {
        fC_res_list[[paste0("Neuron_", i)]] = fC_res_dt
      }
      
    } else if (n_bam == 9) {
      
      if (i <= 3) {
        fC_res_list[[paste0("NES_", i)]] = fC_res_dt
      } else if (7 > i & i > 3) {
        fC_res_list[[paste0("Progenitor_", i)]] = fC_res_dt
      } else if (i > 6) {
        fC_res_list[[paste0("Neuron_", i)]] = fC_res_dt
      }
    }
  }
  

  return(fC_res_list)
  
}

RUN_FOVERLAP <- function(log2_table, gene_DT = merged_fC) {
  
  setkey(log2_table, chr, start, end)
  overlap_DT = foverlaps(gene_DT, log2_table, by.x = c("Chr", "Start", "End"), type="any", which=TRUE)
  overlap_DT = na.omit(overlap_DT)
  overlap_DT = gene_DT[overlap_DT$xid,]
  
  return(overlap_DT)
  
}


LOG2_FILTERING <- function(DT = log2fc_DT, log2_threshold = 1, 
                           bam_indir = "/RNA-BLISS_Data_analysis_Pipeline/Data/RNA_data/RNA_BAM/") {
  
  
  # Filtering for bins with high DSB freq during specific stage of neurodifferetiation
  DT = as.data.table(DT)
  
  # Creating annotation table with gene position of transcript with highest TPM
  # across neurodifferentiation
  
  #bam_indir = "Data/RNA_data/RNA_BAM/"
  bam_files = list.files(bam_indir, pattern = "sortedByCoord.bam")
  
  fC_res_list = RETURN_FEATURE_COUNTS_LIST(keep_End = TRUE, bam_files = bam_files, bam_indir = bam_indir)
  first_res = TRUE
  Overlap_gene_list = list()
  for (i in 1:length(fC_res_list)) {

    if (first_res == TRUE) {
      first_res = FALSE
      merged_fC = fC_res_list[[i]]
      names(merged_fC)[names(merged_fC) == "TPM"] <- "mean_TPM"    
    } else if (first_res == FALSE) {
      fC_res_dt_i = fC_res_list[[i]]
      merged_fC = merge(merged_fC, fC_res_dt_i[,c("GeneID", "TPM")], by = "GeneID")
      merged_fC$mean_TPM = merged_fC$mean_TPM + merged_fC$TPM
      merged_fC = merged_fC[-ncol(merged_fC)]
    }
    
    if (i == 3 | i == 6 | i == 9) {
      
      first_res = TRUE
      
      merged_fC$mean_TPM = merged_fC$mean_TPM / 3
      merged_fC = merged_fC[(merged_fC$mean_TPM >= 1),]
      
      merged_fC$Chr = sub("\\;.*", "", merged_fC$Chr)
      merged_fC$Start = sub("\\;.*", "", merged_fC$Start)
      merged_fC$Start = as.integer(merged_fC$Start)
      merged_fC$End = sub("\\;.*", "", merged_fC$End)
      merged_fC$End = as.integer(merged_fC$End)
      merged_fC$Chr = gsub('chr', '', merged_fC$Chr)
      merged_fC$gene_id = gsub("\\..*","",merged_fC$gene_id)
  
      setDT(merged_fC)
      
      if (i == 3) {
        
        Down_log2 = DT[DT$log2fc_NsPr > log2_threshold &
                         DT$log2fc_PrNe > -log2_threshold & DT$log2fc_PrNe < log2_threshold &
                         DT$log2fc_NsNe > log2_threshold,]
        
        overlap_Down_DT = RUN_FOVERLAP(log2_table = Down_log2, gene_DT = merged_fC)
        Overlap_gene_list[["NES_down"]] = overlap_Down_DT
        
        Up_log2 = DT[DT$log2fc_NsPr < -log2_threshold &
                       DT$log2fc_PrNe > -log2_threshold & DT$log2fc_PrNe < log2_threshold &
                       DT$log2fc_NsNe < -log2_threshold,]
        
        overlap_Up_DT = RUN_FOVERLAP(log2_table = Up_log2, gene_DT = merged_fC)
        Overlap_gene_list[["NES_up"]] = overlap_Up_DT
        
      } else if (i == 6) {
        
        Down_log2 = DT[DT$log2fc_NsPr < -log2_threshold &
                         DT$log2fc_PrNe > log2_threshold &
                         DT$log2fc_NsNe > -log2_threshold & DT$log2fc_NsNe < log2_threshold,]
        
        overlap_Down_DT = RUN_FOVERLAP(log2_table = Down_log2, gene_DT = merged_fC)
        Overlap_gene_list[["Progenitor_down"]] = overlap_Down_DT
        
        Up_log2 = DT[DT$log2fc_NsPr > log2_threshold &
                       DT$log2fc_PrNe < -log2_threshold &
                       DT$log2fc_NsNe > -log2_threshold & DT$log2fc_NsNe < log2_threshold,]
        
        overlap_Up_DT = RUN_FOVERLAP(log2_table = Up_log2, gene_DT = merged_fC)
        Overlap_gene_list[["Progenitor_up"]] = overlap_Up_DT
        
      } else if (i == 9) {
        
        Down_log2 = DT[DT$log2fc_NsPr > -log2_threshold & DT$log2fc_NsPr < log2_threshold &
                         DT$log2fc_PrNe < -log2_threshold &
                         DT$log2fc_NsNe < -log2_threshold,]
        
        overlap_Down_DT = RUN_FOVERLAP(log2_table = Down_log2, gene_DT = merged_fC)
        Overlap_gene_list[["Neuron_down"]] = overlap_Down_DT
        
        Up_log2 = DT[DT$log2fc_NsPr > -log2_threshold & DT$log2fc_NsPr < log2_threshold &
                       DT$log2fc_PrNe > log2_threshold &
                       DT$log2fc_NsNe > log2_threshold,]
        
        overlap_Up_DT = RUN_FOVERLAP(log2_table = Up_log2, gene_DT = merged_fC)
        Overlap_gene_list[["Neuron_up"]] = overlap_Up_DT
        
      }
    }
  }
  
  return(Overlap_gene_list)
  
}

MAP_RNA_BLISS <- function(Overlap_genes, RNA_table_dir, Cross_table_dir) {
  
  venn_list = list()

  for (i in 1:length(Overlap_genes)) {

    if (nrow(Overlap_genes[[i]]) == 0) next
    
    Overlap_suffix = names(Overlap_genes[i])
    Overlap_table = Overlap_genes[[i]]
    RNA_table = paste0(Overlap_suffix, "_dds_raw_counts.rds")
    RNA_table = readRDS(paste0(RNA_table_dir, RNA_table))
                               
    RNA_v = RNA_table$gene_name %in% Overlap_table$gene_id
    Overlap_v = Overlap_table$gene_id %in% RNA_table$gene_name

    RNA_n = length(RNA_v) 
    Overlap_n = length(Overlap_v)
    cross_n = table(RNA_v)["TRUE"]
    
    if (is.na(cross_n) == TRUE) {
      cross_n = 0
    } else if (is.na(cross_n) == FALSE){
      Cross_table = merge(Overlap_table, RNA_table, by.x='gene_id', by.y='gene_name')
      Cross_table = Cross_table[,c(1:5)]
      saveRDS(Cross_table, file = paste0(Cross_table_dir, Overlap_suffix, "_overlapping_genes.rds"))
    }
    
    #grid.newpage()  
    venn.plot <- draw.pairwise.venn(area1      = RNA_n,
                                    area2      = Overlap_n,
                                    cross.area = cross_n,
                                    fill            = c("blue", "red"),
                                    lty             = "solid",
                                    cex             = 2,
                                    cat.cex         = 2,
                                    cat.pos         = c(270, 30),
                                    fontfamily = "sans",
                                    ext.pos         = 30,
                                    ext.dist        = -0.05,
                                    ext.length      = 0.85,
                                    ext.line.lwd    = 2,
                                    ext.line.lty    = "dashed",
                                    ind = FALSE
                                    )
    
    venn_list[[Overlap_suffix]] = venn.plot
    
  }
  
  return(venn_list)
  
}
