# Author: Gustaw Eriksson
# Date: 2020-05-14

# Description: Maps DSB frequency over TSS in genes using RNA-sequencing data.

library("Rsubread")
library("GenomicRanges")
library(GenomicFeatures)
library("data.table")
library("dplyr")
library("ggplot2")
library("cowplot")

## Setting paramters

chr_Y = FALSE
chr_vector = c(1:22, "X")
window_size = 30000
bin_size = 10
local_DSB = FALSE
bin_DSB = TRUE

indir = "/RNA-BLISS_Data_analysis_Pipeline/"
setwd(indir)
dir.create("Output/DSB_mapping", showWarnings = FALSE)
## RNA PART ##
anno = "Data/hg19/gencode.v19.annotation.gtf"

RNA_tables_dir = "Output/RNA_analysis/"
RNA_tables_dir = paste0(RNA_tables_dir, list.files(RNA_tables_dir)[length(list.files(RNA_tables_dir))], "/tables/")

#bam_indir = "/Users/gustaweriksson/Dokument/Bienko-Crosetto/Final_pipeline/Data/RNA_data/RNA_BAM/"
bam_indir = "Data/RNA_data/RNA_BAM/"
bam_files = list.files(bam_indir, pattern = "sortedByCoord.bam")
fC_res_list = list()

## As long as we only have sBLISS for NES, the script is only used for NES and therefore only
## the first three BAM files are used
#for (i in 1:length(bam_files)) {

source("/RNA-BLISS_Data_analysis_Pipeline/Scripts/DSB_mapping/DSB_mapping_tools.R")

#fC_res_list = RETURN_FEATURE_COUNTS_LIST(bam_files = list.files(bam_indir, pattern = "sortedByCoord.bam"), 
#                                         n_bam = 3)

fC_res_list = RETURN_FEATURE_COUNTS_LIST(keep_End = FALSE, bam_files = bam_files, bam_indir = bam_indir, n_bam = 3)

old_script = FALSE
if (old_script == TRUE) {
  for (i in 1:3) {
    
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
    #fC_res_dt = fC_res_dt[(grepl("chrM", fC_res_dt$Chr) == FALSE & grepl("chrY", fC_res_dt$Chr) == FALSE),]
    if (keep_End == FALSE) {
      fC_res_dt = fC_res_dt[c(1:3, 7, 9)]
    } else if (keep_End == TRUE) {
      fC_res_dt = fC_res_dt[c(1:4, 7, 9)]
    }
    
    
    
    fC_res_dt = fC_res_dt[order(fC_res_dt$TPM, decreasing = TRUE),]
    fC_res_dt = fC_res_dt[(!duplicated(fC_res_dt$gene_id) == TRUE),]
    fC_res_dt = fC_res_dt[fC_res_dt$TPM >= 1,]
    
    if (length(bam_files) == 6 | length(bam_files) == 3) {
      
      if (i <= 3) {
        fC_res_list[[paste0("NES_", i)]] = fC_res_dt
      } else if (i > 3) {
        fC_res_list[[paste0("Neuron_", i)]] = fC_res_dt
      }
      
    } else if (length(bam_files) == 9) {
      
      if (i <= 3) {
        fC_res_list[[paste0("NES_", i)]] = fC_res_dt
      } else if (7 > i & i > 3) {
        fC_res_list[[paste0("Progenitor_", i)]] = fC_res_dt
      } else if (i > 6) {
        fC_res_list[[paste0("Neuron_", i)]] = fC_res_dt
      }
    }
  }
}

first_res = TRUE
merged_fC_res_list = list()
top_5_list = list()
median_list = list()
bottom_5_list = list()
Up_exp_list = list()
Down_exp_list = list()

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
    
    suffix_i = sub("\\_.*", "", names(fC_res_list[i]))
    merged_fC$Chr = sub("\\;.*", "", merged_fC$Chr)
    merged_fC$Start = sub("\\;.*", "", merged_fC$Start)
    merged_fC$Chr = gsub('chr', '', merged_fC$Chr)
    merged_fC$gene_id = gsub("\\..*","",merged_fC$gene_id)
    names(merged_fC)[names(merged_fC) == "Start"] <- "Pos"
    
    ### FILTER OUT TOP AND BOTTOM 5% DEPENDING ON TPM ###
    merged_fC_TPM = merged_fC
    TPM_5 = quantile(merged_fC$mean_TPM, c(0.05, 0.475, 0.525, 0.95))
    
    bottom_5_DT = subset(merged_fC, mean_TPM <= TPM_5[[1]])
    bottom_5_DT = bottom_5_DT[-c(1,4,5)]
    median_DT = subset(merged_fC, mean_TPM >= TPM_5[[2]] & mean_TPM <= TPM_5[[3]])
    median_DT = median_DT[-c(1,4,5)]
    top_5_DT = subset(merged_fC, mean_TPM >= TPM_5[[4]])
    top_5_DT = top_5_DT[-c(1,4,5)]
    
    if (i == 3) {
      bottom_5_list[["NES"]] = bottom_5_DT
      median_list[["NES"]] = median_DT
      top_5_list[["NES"]] = top_5_DT
      
      RNA_table = readRDS(paste0(RNA_tables_dir, list.files(RNA_tables_dir, pattern = "Up_NES_dds_raw_counts.rds")))
      NES_up_DT = subset(merged_fC, gene_id %in% RNA_table$gene_name)
      NES_up_DT = NES_up_DT[-c(1,4,5)]
      Up_exp_list[["NES"]] = NES_up_DT
      
      RNA_table = readRDS(paste0(RNA_tables_dir, list.files(RNA_tables_dir, pattern = "Down_NES_dds_raw_counts.rds")))
      NES_down_DT = subset(merged_fC, gene_id %in% RNA_table$gene_name)
      NES_down_DT = NES_down_DT[-c(1,4,5)]
      Down_exp_list[["NES"]] = NES_down_DT
      
      merged_fC = merged_fC[-c(1,4,5)]
      merged_fC_res_list[["NES"]] = merged_fC
      
    } else if (i == 6) {
      bottom_5_list[["Progenitor"]] = bottom_5_DT
      median_list[["Progenitor"]] = median_DT
      top_5_list[["Progenitor"]] = top_5_DT
      
      RNA_table = readRDS(paste0(RNA_tables_dir, list.files(RNA_tables_dir, pattern = "Up_Progenitor_dds_raw_counts.rds")))
      Progenitor_up_DT = subset(merged_fC, gene_id %in% RNA_table$gene_name)
      Progenitor_up_DT = Progenitor_up_DT[-c(1,4,5)]
      Up_exp_list[["Progenitor"]] = Progenitor_up_DT
      
      RNA_table = readRDS(paste0(RNA_tables_dir, list.files(RNA_tables_dir, pattern = "Down_Progenitor_dds_raw_counts.rds")))
      Progenitor_down_DT = subset(merged_fC, gene_id %in% RNA_table$gene_name)
      Progenitor_down_DT = Progenitor_down_DT[-c(1,4,5)]
      Down_exp_list[["Progenitor"]] = Progenitor_down_DT
      
      merged_fC = merged_fC[-c(1,4,5)]
      merged_fC_res_list[["Progenitor"]] = merged_fC
    } else if (i == 9) {
      bottom_5_list[["Neuron"]] = bottom_5_DT
      median_list[["Neuron"]] = median_DT
      top_5_list[["Neuron"]] = top_5_DT
      
      RNA_table = readRDS(paste0(RNA_tables_dir, list.files(RNA_tables_dir, pattern = "Up_Neuron_dds_raw_counts.rds")))
      Neuron_up_DT = subset(merged_fC, gene_id %in% RNA_table$gene_name)
      Neuron_up_DT = Neuron_up_DT[-c(1,4,5)]
      Up_exp_list[["Neuron"]] = Neuron_up_DT
      
      RNA_table = readRDS(paste0(RNA_tables_dir, list.files(RNA_tables_dir, pattern = "Down_Neuron_dds_raw_counts.rds")))
      Neuron_down_DT = subset(merged_fC, gene_id %in% RNA_table$gene_name)
      Neuron_down_DT = Neuron_down_DT[-c(1,4,5)]
      Down_exp_list[["Neuron"]] = Neuron_down_DT
      
      merged_fC = merged_fC[-c(1,4,5)]
      merged_fC_res_list[["Neuron"]] = merged_fC
    }
  }
}

## BLISS PART ##

source("/RNA-BLISS_Data_analysis_Pipeline/Scripts/Data_handling/core_functions.R")

## Load the BLISS bed files, expand when more BLISS data is available
bed_indir = "/RNA-BLISS_Data_analysis_Pipeline//Data/BLISS_bed/sBLISS/"
load_DSB_DT_list = load_bed_from_dir(indir = bed_indir)

COUNT_DSB_TSS <- function(Position_DT, DSB_DT_list = load_DSB_DT_list) {
  
  # Generating DT for DSB pos up and downstream from TSS
  count_DT = data.table(Pos = as.character(c((-window_size/2):-1, 0, 1:(window_size/2))), 
                        DSB = as.numeric(rep(0, window_size + 1)))
  
  for (i in 1:length(DSB_DT_list)) {
    
    first_DSB_DT = TRUE
    
    for (j in 1:length(DSB_DT_list[[i]])) {

      print(paste0("FILE NR ", j, "/", length(DSB_DT_list[[i]])))
      DSB_DT = as.data.table(DSB_DT_list[[i]][[j]])
      names(DSB_DT)[names(DSB_DT) == "seqnames"] <- "Chr"
      names(DSB_DT)[names(DSB_DT) == "start"] <- "Pos"
      DSB_DT$Chr = as.character(DSB_DT$Chr)
      
      DSB_DT = DSB_DT[grepl("chr", DSB_DT$Chr)]
      if (chr_Y == FALSE) {
        DSB_DT = DSB_DT[!grepl("chrY", DSB_DT$Chr)]
      }
      DSB_DT$Chr = gsub('chr', '', DSB_DT$Chr)
      DSB_DT = DSB_DT[,-c(3:5)]
      
      ## Loop over each chr
      for (chr_n in chr_vector) {
        
        print(chr_n)
        
        ## Create vector with TSS position from the RNA data for the current chr
        TSS_pos_vector = Position_DT[which(Position_DT$Chr == chr_n), ]
        TSS_pos_vector$Pos = as.numeric(TSS_pos_vector$Pos)
        TSS_pos_vector = TSS_pos_vector[order(TSS_pos_vector$Pos),]
        
        ## Grep BLISS with current chr
        BLISS_DT = DSB_DT[which(DSB_DT$Chr == chr_n), ]
        
        for (pos in TSS_pos_vector$Pos) {
          
          min_pos = pos - (window_size/2 + bin_size/2)
          max_pos = pos + (window_size/2 + bin_size/2)
          
          DSB_pos = BLISS_DT[which(BLISS_DT$Pos >= min_pos & BLISS_DT$Pos <= max_pos)]
          DSB_pos = DSB_pos[,-1]
          
          count_DT$DSB[(window_size/2)+1] = count_DT$DSB[(window_size/2) + 1] + 
            sum(DSB_pos$count[which(DSB_pos$Pos >= pos - (bin_size/2) & DSB_pos$Pos <= pos + (bin_size/2))])
          
          
          TSS_count = count_DT$DSB[(window_size/2) + 1] + 
            sum(DSB_pos$count[which(DSB_pos$Pos >= pos - (bin_size/2) & DSB_pos$Pos <= pos + (bin_size/2))])

          DSB_pos = DSB_pos[!(DSB_pos$Pos >= pos - (bin_size/2) & DSB_pos$Pos <= pos + (bin_size/2))]
          
          DSB_pos$Pos = DSB_pos$Pos - pos
          DSB_pos$Pos[which(DSB_pos$Pos < 0)] = DSB_pos$Pos[which(DSB_pos$Pos < 0)] + bin_size/2
          DSB_pos$Pos[which(DSB_pos$Pos > 0)] = DSB_pos$Pos[which(DSB_pos$Pos > 0)] - bin_size/2
          DSB_pos$Pos = as.character(DSB_pos$Pos)
          
          count_DT = left_join(count_DT, DSB_pos, by = "Pos")
          count_DT[is.na(count_DT)] = 0
          count_DT$DSB = count_DT$DSB + count_DT$count
          count_DT = count_DT[,-3]
          
        }
        
      }
      
    }
    
  }
  
  n_breaks = window_size/bin_size
  count_DT$Pos = c(rep((-n_breaks/2):-1, each = bin_size), 0, rep(1:(n_breaks/2), each = bin_size))
  count_DT = count_DT %>% group_by(Pos) %>% mutate(., Normalised_DSB = mean(DSB))
  count_DT$Normalised_DSB[(window_size/2)+1] = count_DT$Normalised_DSB[(window_size/2)+1]/bin_size
  count_DT = count_DT[!duplicated(count_DT$Pos),]
  #count_DT$Normalised_DSB = count_DT$Normalised_DSB / bin_size
  count_DT$Pos = as.numeric(count_DT$Pos)
  count_DT$Pos = as.factor(count_DT$Pos)
  #count_DT_spline = as.data.frame(spline(count_DT$Pos, count_DT$Normalised_DSB))
  
  return(count_DT)
  
}

PLOT_COUNT <- function(count_DT, count_DT_2 = NULL, count_DT_3 = NULL, window_n) {
  
  n_breaks = window_size/bin_size
  window_n = window_n/1000
  window_label = c(paste0(-window_n/2, "Kb"), paste0(-window_n/4, "Kb"),
                   "TSS", paste0(window_n/4, "Kb"), paste0(window_n/2, "Kb"))
  count_DT_spline = as.data.frame(spline(count_DT$Pos, count_DT$Normalised_DSB))
  
  if (is.null(count_DT_2) == TRUE) {
    
    plot_TSS_rough = ggplot(count_DT, aes(x = Pos, y = Normalised_DSB, group = 1)) +
      geom_line(data = count_DT_spline, aes(x = x, y = y), size = 1) + 
      theme(plot.title = element_text(size = 15), axis.text.x = element_text(angle = 45, hjust = 1)) +
      ylab("Normalised DSB count") + theme_cowplot() + ggtitle("DSB in proximity of TSS \nBin size = 10bp") +
      geom_vline(xintercept = n_breaks/2 + 1, colour = "red", linetype = "dashed", size = 0.5) +
      scale_x_discrete(name = "Relative position", breaks = c(-n_breaks/2, -n_breaks/4, 0, n_breaks/4 ,n_breaks/2),
                       labels = window_label)

    ggsave(paste0("Output/DSB_mapping/All_", window_n, "Kb_rough_DSB_TSS_plot.pdf"), plot = plot_TSS_rough,
           device = "pdf")
    
    plot_TSS = ggplot(count_DT, aes(x = Pos, group = 1)) +
      geom_smooth(data = count_DT, aes(y = Normalised_DSB), method = "loess", colour = "black", se = FALSE) +
      ylab("Normalised DSB count") + theme_cowplot() + ggtitle("DSB in proximity of TSS \nBin size = 10bp") +
      scale_x_discrete(name = "Relative position", breaks = c(-n_breaks/2, -n_breaks/4, 0, n_breaks/4 ,n_breaks/2), 
                       labels = window_label) +
      theme(plot.title = element_text(size = 15), axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_vline(xintercept = n_breaks/2 + 1, colour = "black", linetype = "dashed", size = 0.5)
    
    ggsave(paste0("Output/DSB_mapping/All_", window_n, "Kb_smooth_DSB_TSS_plot.pdf"), plot = plot_TSS,
           device = "pdf")
    
    plot_list = list()
    plot_list[["rough"]] = plot_TSS_rough
    plot_list[["smooth"]] = plot_TSS
    
    
  } else if (is.null(count_DT_2) == FALSE & is.null(count_DT_3) == TRUE) {
    
    count_DT_2_spline = as.data.frame(spline(count_DT_2$Pos, count_DT_2$Normalised_DSB))
    count_DT_spline$y_2 = count_DT_2_spline$y
    
    plot_TSS = ggplot(count_DT, aes(x = Pos, group = 1)) +
      geom_smooth(data = count_DT, aes(y = Normalised_DSB), method = "loess", colour = "red") +
      geom_smooth(data = count_DT_2, aes(y = Normalised_DSB), method = "loess", colour = "blue") +
      ylab("Normalised DSB count") + theme_cowplot() + ggtitle("DSB in proximity of TSS \nBin size = 10bp") +
      scale_x_discrete(name = "Relative position", breaks = c(-n_breaks/2, -n_breaks/4, 0, n_breaks/4 ,n_breaks/2), 
                       labels = window_label) +
      theme(plot.title = element_text(size = 15), axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_vline(xintercept = n_breaks/2 + 1, colour = "black", linetype = "dashed", size = 0.5)
    
    ggsave(paste0("Output/DSB_mapping/2_levels_", window_n, "_smooth_DSB_TSS_plot.pdf"), plot = plot_TSS,
           device = "pdf")
    
    return(plot_TSS)
    
  } else if (is.null(count_DT_2) == FALSE & is.null(count_DT_3) == FALSE) {
    
    plot_TSS = ggplot(count_DT, aes(x = Pos, group = 1)) +
      geom_smooth(data = count_DT, aes(y = Normalised_DSB), method = "loess", colour = "red", se = FALSE) +
      geom_smooth(data = count_DT_2, aes(y = Normalised_DSB), method = "loess", colour = "blue", se = FALSE) +
      geom_smooth(data = count_DT_3, aes(y = Normalised_DSB), method = "loess", colour = "green", se = FALSE) +
      ylab("Normalised DSB count") + theme_cowplot() + ggtitle("DSB in proximity of TSS \nBin size = 10bp") +
      scale_x_discrete(name = "Relative position", breaks = c(-n_breaks/2, -n_breaks/4, 0, n_breaks/4 ,n_breaks/2), 
                       labels = window_label) +
      theme(plot.title = element_text(size = 15), axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_vline(xintercept = n_breaks/2 + 1, colour = "black", linetype = "dashed", size = 0.5)

    ggsave(paste0("Output/DSB_mapping/3_levels_", window_n, "_smooth_DSB_TSS_plot.pdf"), plot = plot_TSS,
           device = "pdf")
    
    return(plot_TSS)
    
  }
  
}

Main_count = COUNT_DSB_TSS(Position_DT = merged_fC_res_list[[1]])
top_5_count = COUNT_DSB_TSS(Position_DT = top_5_DT)
median_count = COUNT_DSB_TSS(Position_DT = median_DT)
bottom_5_count = COUNT_DSB_TSS(Position_DT = bottom_5_DT)

window_suffix = window_size/1000

saveRDS(Main_count, file = paste0("Output/DSB_mapping/Main_", window_suffix, "Kb_count.rds"))
saveRDS(top_5_count, file = paste0("Output/DSB_mapping/top_5_", window_suffix, "Kb_count.rds"))
saveRDS(median_count, file = paste0("Output/DSB_mapping/median_", window_suffix, "Kb_count.rds"))
saveRDS(bottom_5_count, file = paste0("Output/DSB_mapping/bottom_5_", window_suffix, "Kb_count.rds"))

plot_main_list = PLOT_COUNT(count_DT = Main_count, window_n = window_size)
plot_5_median = PLOT_COUNT(count_DT = top_5_count, count_DT_2 = bottom_5_count, count_DT_3 = median_count, window_n = window_size)

