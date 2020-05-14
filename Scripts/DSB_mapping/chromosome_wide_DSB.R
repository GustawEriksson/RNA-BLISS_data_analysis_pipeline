# Author: Gustaw Eriksson
# Data: 18-12-2019
# Description: Mapping DSBs across the chromosomes and calculating log2fc of bins between groups

run_chr_dsb_BLISS = chr_WIDE_DSB(indir = "/RNA-BLISS_Data_analysis_Pipeline/Output/BLISS", 
                                 cell_types = c("NES", "Progenitor", "Neuron"), BLISS_run = "B138", window_size = "150Kb",
                                 merge_types = TRUE,
                                 log2fc_analysis = TRUE, global_DSB = FALSE)

run_chr_dsb_sBLISS = chr_WIDE_DSB(indir = "/RNA-BLISS_Data_analysis_Pipeline/Output/sBLISS/", 
                                  cell_types = "NES", BLISS_run = "B194", window_size = "10Kb",
                                  merge_types = FALSE, 
                                  log2fc_analysis = TRUE)

chr_WIDE_DSB <- function(indir = "/RNA-BLISS_Data_analysis_Pipeline/Output/", 
                         cell_types = c("NES", "Progenitor", "Neuron"), BLISS_run, window_size,
                         merge_types = FALSE, 
                         log2fc_analysis = TRUE, 
                         log2_limit = log2(1.5),
                         global_DSB = FALSE) {
  
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
  
  source("/RNA-BLISS_Data_analysis_Pipeline/Scripts/DSB_mapping/DSB_mapping_tools.R")
  
  #source("/Users/gustaweriksson/Dokument/Bienko-Crosetto/Final_pipeline/Scripts/create_merged_bedtools_tables/bedtools_create_merged_count_table.R")
  setwd(indir)
  dir.create("BLISS", showWarnings = FALSE)
  indir = paste0(indir, "BLISS/")
  setwd(indir)
  dir.create(BLISS_run, showWarnings = FALSE)
  dir.create(paste0(BLISS_run, "/", window_size), showWarnings = FALSE)
  dir.create(paste0(BLISS_run, "/", window_size, "/chr_plots"), showWarnings = FALSE)
  dir.create(paste0(BLISS_run, "/", window_size, "/Log2FC_analysis"), showWarnings = FALSE)
  chr_vector = c(1:22, "X")
  cell_line_plot_list = list()
  cell_normalised_line_plot_list = list()
  first_merge = TRUE
  
  for (type in cell_types){
    
    dir.create(paste0(BLISS_run, "/", window_size, "/", type, "_chr_plots"), showWarnings = FALSE)
    # Load in the raw count or normalised count table
    load_indir = paste0(indir, BLISS_run, "/", window_size, "/Count_tables/")
    
    load_file = paste0(load_indir, list.files(path = load_indir, pattern = paste0(type,"_", window_size, "_bins_raw_count.csv")))
    raw_DT = EDIT_COUNT_DT(count_dt_dir = load_file)
    if (first_merge == TRUE){
      first_merge = FALSE
      sum_count = sum(raw_DT$Sum_Count)
      merge_raw = raw_DT[-c(4:(length(raw_DT)-1))] %>% mutate(Sum_Count = (Sum_Count/sum_count)*1000000)
      
      plotting_DSB_chr = PLOT_DSB_CPM_ALL_chr(CPM_DT = merge_raw, plot_title = type, BLISS_run = BLISS_run, window_size = window_size)
      
      names(merge_raw)[length(merge_raw)] = paste0(type, "_Normalised_Sum_Count")

    } else if (merge_types == TRUE & first_merge == FALSE) {
      sum_count = sum(raw_DT$Sum_Count)
      norm_raw_DT = raw_DT %>% mutate(Sum_Count = (Sum_Count/sum_count)*1000000)
      plotting_DSB_chr = PLOT_DSB_CPM_ALL_chr(CPM_DT = norm_raw_DT, plot_title = type, BLISS_run = BLISS_run, window_size = window_size)
      merge_raw$tmp = norm_raw_DT$Sum_Count
      names(merge_raw)[length(merge_raw)] = paste0(type, "_Normalised_Sum_Count")
    }
    
    # Subset DT depending on chr number
    plot_line_single_chr_list = list()
    plot_bar_single_chr_list = list()
    
    for (chr_i in chr_vector){
      chr_count_DT = merge_raw[which(merge_raw$chr == chr_i), ]
      
      names(chr_count_DT)[ncol(chr_count_DT)] = "CPM"
      
      min_end_int = min(chr_count_DT$end)
      mid_end_int = chr_count_DT[round(nrow(chr_count_DT)/2) - 1, 3]
      mid_end_Mb = round(mid_end_int/1000000)
      max_end_int = max(chr_count_DT$end)
      max_end_Mb = round(max_end_int/1000000)
      
      chr_count_DT$end <- as.factor(chr_count_DT$end)
      
      chr_line_plot = ggplot(data = chr_count_DT, aes(x = end, y = CPM, group = 1)) +
        geom_line() + geom_hline(yintercept = median(chr_count_DT$CPM), color="red") +
        theme_cowplot() + xlab("") + ylab("") + ggtitle(paste0("chr ", chr_i)) + ylim(0, 200) +
        theme(plot.title = element_text(hjust = 0.5, size = 10)) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
        theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
      
      chr_line_plot_single = ggplot(data = chr_count_DT, aes(x = end, y = CPM, group = 1)) +
        geom_area(color = 'black') + scale_x_discrete(name="Chromosome position (Mb)", 
                                                      breaks=c(min_end_int, mid_end_int, max_end_int), 
                                                      labels=c("0", mid_end_Mb, max_end_Mb)) +
        geom_hline(yintercept = median(chr_count_DT$CPM), color="red") +
        ylab("DSB counts per million") + ggtitle(paste0("chr ", chr_i, ", bin size = ", window_size)) + 
        ylim(0, max(chr_count_DT$CPM)*1.1) + theme_cowplot() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

      ggsave(filename = paste0("chr_",chr_i, "_", window_size, "_", type, "_", BLISS_run, "_chr_line_plot.pdf"), plot = chr_line_plot_single, device = "pdf", 
             dpi = "retina", scale = 1, limitsize = FALSE, path = paste0(BLISS_run, "/", window_size, "/", type, "_chr_plots"))
        
      plot_line_single_chr_list[[chr_i]] = chr_line_plot
      
      chr_bar_plot = ggplot(data = chr_count_DT, aes(x = end, y = Sum_Count)) +
        geom_bar(stat="identity", color = "black") + geom_hline(yintercept = median(chr_count_DT$Sum_Count), color="red") + scale_x_continuous(labels = comma) +
        theme_cowplot() + xlab("chromosome bp") + ylab("Mean bin DSB count") + ggtitle(paste0("Mean DSB count per bin across chromosome ", chr_i))
    }

      cell_line_plot_list[[type]] = plot_line_single_chr_list
    
    
  }
  
  cols = 4
  rows = 6
  
    
  if (is.element("NES", cell_types) == TRUE) {
    NES_arrange = ggarrange(plotlist = cell_line_plot_list$NES, ncol=cols, nrow = rows) + ggtitle("NES")
    NES_arrange = annotate_figure(NES_arrange, top = paste0("DSB frequency across NES, ", BLISS_run,". \n Bin size = ", window_size))
    ggsave(filename = paste0(window_size, "_", BLISS_run, "_NES_chr_line_plot.pdf"), plot = NES_arrange, device = "pdf", 
           dpi = "retina", scale = 4, limitsize = FALSE, path = paste0(BLISS_run, "/", window_size, "/NES_chr_plots"))
  }
  
  if (is.element("Progenitor", cell_types) == TRUE) {
    Pr_arrange = ggarrange(plotlist = cell_line_plot_list$Progenitor, ncol=cols, nrow = rows)
    Pr_arrange = annotate_figure(Pr_arrange, top = paste0("DSB frequency across progenitors, ", BLISS_run,". \n Bin size = ", window_size))
    ggsave(filename = paste0(window_size, "_", BLISS_run, "_Progenitor_chr_line_plot.pdf"), plot = Pr_arrange, device = "pdf", 
           dpi = "retina", scale = 4, limitsize = FALSE, path = paste0(BLISS_run, "/", window_size, "/Progenitor_chr_plots"))
  }
  
  if (is.element("Neuron", cell_types) == TRUE) {
    Ne_arrange = ggarrange(plotlist = cell_line_plot_list$Neuron, ncol=cols, nrow = rows)
    Ne_arrange = annotate_figure(Ne_arrange, top = paste0("DSB frequency across neurons, ", BLISS_run,". \n Bin size = ", window_size))
    ggsave(filename = paste0(window_size, "_", BLISS_run, "_Neuron_chr_line_plot.pdf"), plot = Ne_arrange, device = "pdf", 
           dpi = "retina", scale = 4, limitsize = FALSE, path = paste0(BLISS_run, "/", window_size, "/Neuron_chr_plots"))
  }
  
  if (merge_types == TRUE) {
    plot_melt_line_single_chr_list = list()
    plot_melt_bar_single_chr_list = list()

    for (chr_i in chr_vector){
      chr_merge_count_DT = merge_raw[which(merge_raw$chr == chr_i), ]
      melt_chr_merge_count_DT = chr_merge_count_DT[,-c(1:2)]
      names(melt_chr_merge_count_DT) = c("end", cell_types)
      melt_chr_merge_count_DT$end = melt_chr_merge_count_DT$end/1000000
      
      melt_chr_merge_count_DT = melt(melt_chr_merge_count_DT, id.vars = 'end', variable.name = 'Cell_type')
      
      melt_line_plot = ggplot(melt_chr_merge_count_DT, aes(x = end, y = value, colour = Cell_type)) + 
        geom_line(aes(colour = Cell_type), size = 0.8) + theme_cowplot() + 
        theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), 
              axis.ticks.y=element_blank()) + xlab("") + ylab("") + ggtitle(paste0("chr ", chr_i)) +
        theme(plot.title = element_text(hjust = 0.5, size = 30)) + 
        theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) +
        guides(colour = guide_legend(override.aes = list(size=10)))
      
      melt_line_plot_single = ggplot(melt_chr_merge_count_DT, aes(x = end, y = value, colour = Cell_type)) + 
        geom_line(aes(colour = Cell_type), size = 1) + theme_cowplot() +
        xlab("Chromosome position (Mb)") + ylab("DSB count per million") + ggtitle(paste0("chr ", chr_i)) +
        theme(plot.title = element_text(hjust = 0.5, size = 20))
      
      ggsave(filename = paste0(window_size, "_", BLISS_run, "_chr_", chr_i, "_line_plot.pdf"), plot = melt_line_plot_single, device = "pdf", 
             dpi = "retina", scale = 1, limitsize = FALSE, path = paste0(BLISS_run, "/", window_size, "/chr_plots"))
      
      melt_bar_plot = ggplot(melt_chr_merge_count_DT, aes(x = end, y = value, colour = Cell_type)) + 
        geom_bar(position = "dodge", stat = "identity") + theme_cowplot() + 
        theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), 
              axis.ticks.y=element_blank()) + xlab("") + ylab("") + ggtitle(paste0("chr ", chr_i)) +
        theme(plot.title = element_text(hjust = 0.5, size = 30)) + 
        theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20)) +
        guides(colour = guide_legend(override.aes = list(size=10)))
      
      plot_melt_line_single_chr_list[[chr_i]] = melt_line_plot
      plot_melt_bar_single_chr_list[[chr_i]] = melt_bar_plot
      
    }
    
    All_line_arrange = ggarrange(plotlist = plot_melt_line_single_chr_list, ncol=cols, nrow = rows,
                                 legend = "top", common.legend = TRUE) 
    All_line_arrange = annotate_figure(All_line_arrange, top = text_grob(paste0("All cell types, ", window_size), size = 40))
    ggsave(filename = paste0(window_size, "_All_Cell_chr_line_plot.pdf"), plot = All_line_arrange, device = "pdf", 
           dpi = "retina", scale = 4, limitsize = FALSE, path = paste0(BLISS_run, "/", window_size, "/chr_plots"))
    All_bar_arrange = ggarrange(plotlist = plot_melt_bar_single_chr_list, ncol=cols, nrow = rows,
                                legend = "top", common.legend = TRUE)
    All_bar_arrange = annotate_figure(All_bar_arrange, top = text_grob(paste0("All cell types, ", window_size), size = 40))
    ggsave(filename = paste0(window_size, "_All_Cell_chr_bar_plot.pdf"), plot = All_bar_arrange, device = "pdf", 
           dpi = "retina", scale = 4, limitsize = FALSE, path = paste0(BLISS_run, "/", window_size, "/chr_plots"))
  }
  
  if (log2fc_analysis == TRUE) {
    
    log2fc_vector = c("Log2FC NsNe", "Log2FC NsPr", "Log2FC PrNe")
    
    log2fc_DT = merge_raw
    log2fc_DT$log2fc_NsNe = log2((merge_raw$Neuron_Normalised_Sum_Count+1)/(merge_raw$NES_Normalised_Sum_Count+1))
    log2fc_DT$log2fc_NsPr = log2((merge_raw$Progenitor_Normalised_Sum_Count+1)/(merge_raw$NES_Normalised_Sum_Count+1))
    log2fc_DT$log2fc_PrNe = log2((merge_raw$Neuron_Normalised_Sum_Count+1)/(merge_raw$Progenitor_Normalised_Sum_Count+1))
    log2fc_DT = log2fc_DT[,-c(4:6)]
    
    saveRDS(log2fc_DT, file = paste0(BLISS_run, "/", window_size, "/Log2FC_analysis/Log2FC_DT.rds"))
    
    Overlap_gene_list = LOG2_FILTERING(DT = log2fc_DT, 
                                       log2_threshold = log2_limit, 
                                       bam_indir = "/RNA-BLISS_Data_analysis_Pipeline/Data/RNA_data/RNA_BAM/")
    
    map_RNA_BLISS_list = MAP_RNA_BLISS(Overlap_genes = Overlap_gene_list, 
                                       RNA_table_dir = "/RNA-BLISS_Data_analysis_Pipeline/Output/RNA_analysis/LFC_1_TPM_1/tables/",
                                       Cross_table_dir = paste0(indir, BLISS_run, "/", window_size, "/Log2FC_analysis/"))
    
    for (i in 1:length(map_RNA_BLISS_list)) {
      
      pdf(paste0(BLISS_run, "/", window_size, "/Log2FC_analysis/", names(map_RNA_BLISS_list[i]), "_venn_diagram.pdf"))
      
      grid.draw(map_RNA_BLISS_list[[i]])
      
      dev.off()
      
    }
    
    #Checking_bin_genes = CALLING_BIN_GENE(selected_bins = log2fc_DT, log2FC_threshold = log2(1.5))

    for (log2fc_i in 1:length(log2fc_vector)) {
      
      current_log2fc_type = log2fc_DT[,c(1:3, (log2fc_i+3))]
      names(current_log2fc_type)[4] = "Log2fc"
      
      if (global_DSB == TRUE) {
        
        # Make all log2fc values absolute and thereafter select which are outliers of the 0.95 limit, single tail.
        global_DT = current_log2fc_type
        global_DT$end = seq(from = min(merge_raw$end), to = nrow(merge_raw) * min(merge_raw$end), by = min(merge_raw$end))
        
        global_DT$end = as.factor(global_DT$end)
        global_DT$Log2fc = as.numeric(global_DT$Log2fc)
        
        global_DT_chr_line_plot = ggplot(data = global_DT, aes(x = end, y = Log2fc, group = 1)) +
          geom_area() + ylab(paste0(log2fc_vector[[log2fc_i]], " CPM of DSB")) + 
          geom_smooth(method = loess, se = FALSE) +
          ggtitle(paste0(log2fc_vector[[log2fc_i]], " of CPM of DSB in chr ", chr_vector[[log2fc_j]], "\nBin size = ", window_size)) + 
          xlab("Chromosome position") + theme(plot.title = element_text(size = 20))+
          geom_hline(yintercept = c(0, log2(log2_limit), -log2(log2_limit)), color=c("black", "red", "red"), linetype = c("solid", "dashed", "dashed")) + 
          scale_x_discrete(name = "Chromosome position", breaks = c(min_end_int, mid_end_int, max_end_int), 
                           labels = c("0Mb", paste0(mid_end_Mb, "Mb"), paste0(max_end_Mb, "Mb"))) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5)) + theme_cowplot()
        
        ggsave(filename = paste0(log2fc_vector[[log2fc_i]], "_Global_plot.pdf"), plot = global_DT_chr_line_plot,  device = "pdf", 
               dpi = "retina", limitsize = FALSE, path = paste0(BLISS_run, "/", window_size, "/Log2FC_analysis"))
        
      }
      
      for (log2fc_j in 1:length(chr_vector)) {
        
        current_log2fc = current_log2fc_type[which(current_log2fc_type$chr == chr_vector[[log2fc_j]]),]
        
        current_log2fc = current_log2fc[,-c(1,2)]
        current_log2fc$end = (current_log2fc$end)/1000000

        min_end_int = min(current_log2fc$end)
        mid_end_int = current_log2fc[nrow(current_log2fc)/2,1]
        max_end_int = max(current_log2fc$end)

        current_log2fc$end = as.factor(current_log2fc$end)
        current_log2fc$Log2fc = as.numeric(current_log2fc$Log2fc)
        
        log2fc_chr_line_plot = ggplot(data = current_log2fc, aes(x = end, y = Log2fc, group = 1)) +
          geom_area() + theme_cowplot() +
          ggtitle(paste0(log2fc_vector[[log2fc_i]], " of CPM of DSB in chr ", chr_vector[[log2fc_j]], "\nBin size = ", window_size)) +
          xlab("Chromosome position") + ylab("Log2 fold change") +
          theme(plot.title = element_text(size = 15), axis.text.x = element_text(angle = 45, hjust = 1)) +
          ylim(min(current_log2fc[2]) * 1.01, max(current_log2fc[2]) * 1.01) +
          geom_hline(yintercept = c(0, log2(log2_limit), -log2(log2_limit)), color=c("black", "red", "red"), linetype = c("solid", "dashed", "dashed")) +
          scale_x_discrete(name = "Chromosome position (Mb)", breaks = c(min_end_int, mid_end_int, max_end_int), 
                           labels = c("0", mid_end_int, round(max_end_int)))
        
        ggsave(filename = paste0(log2fc_vector[[log2fc_i]], "_chr_", chr_vector[[log2fc_j]], "_plot.pdf"), plot = log2fc_chr_line_plot,  device = "pdf", 
               dpi = "retina", limitsize = FALSE, path = paste0(BLISS_run, "/", window_size, "/Log2FC_analysis"))
  
        
        
      }
      
      current_log2fc_type$chr <- factor(current_log2fc_type$chr, levels=unique(current_log2fc_type$chr))
      
      stat1 <- current_log2fc_type %>%
        group_by(chr) %>%
        summarise(mean = mean(Log2fc), int1 = as.vector(CI(Log2fc))[1], int2 = as.vector(CI(Log2fc))[3]) 
      
      #To check how reliable Your intervals will be (rather unimportant with sequencing current_log2fc_type (huge n)) - are se values 2 times lower than IC? - what is the confidence coefficient?
      
      #Calculate SE
      stat2 <- current_log2fc_type %>%
        group_by(chr) %>%
        summarise(mean = mean(Log2fc), int1 = as.vector(STDERR(Log2fc))[1], int2 = as.vector(STDERR(Log2fc))[3]) 
      
      #Boxplot
      ggplot(current_log2fc_type) +
        geom_boxplot(aes(y = Log2fc, x = chr)) +
        geom_hline(yintercept = c(0, log2(log2_limit), -log2(log2_limit)), color=c("black", "red", "red"), linetype = c("solid", "dashed", "dashed")) +
        geom_errorbar(data = stat1, aes(x = chr, ymin = int1 - mean, ymax = int2 - mean), color = "red") +
        stat_summary(fun.y = "mean", geom="point", aes(y = Log2fc, x = chr)) +
        scale_x_discrete(limits = c(levels(current_log2fc_type$chr))) +
        ggtitle(paste0(log2fc_vector[[log2fc_i]], " of CPM of DSB all across chromosomes \nBin size = ", window_size)) +
        xlab("Chromosome") + ylab("Log2 fold change") + theme(plot.title = element_text(size = 15)) + theme_cowplot()
      
      ggsave(paste0(log2fc_vector[log2fc_i], "_All_chr_barplot.pdf"), width = 10, height = 7,  dpi = 500, path = paste0(BLISS_run, "/", window_size, "/Log2FC_analysis"))
      
      #Take a closer look at the distribution with violin plots!
      ggplot(current_log2fc_type) +
        geom_violin(aes(y = Log2fc, x = chr)) +
        geom_hline(yintercept = c(0, log2(log2_limit), -log2(log2_limit)), color=c("black", "red", "red"), linetype = c("solid", "dashed", "dashed")) +
        geom_errorbar(data = stat1, aes(x = chr, ymin = int1 - mean, ymax = int2 - mean), color = "red") +
        stat_summary(fun.y = "mean", geom="point", aes(y = Log2fc, x = chr), color = "red") +
        stat_summary(fun.y = "median", geom="point", aes(y = Log2fc, x = chr), color = "blue") +
        scale_x_discrete(limits = c(levels(current_log2fc_type$chr))) +
        ggtitle(paste0(log2fc_vector[[log2fc_i]], " of CPM of DSB across all chromosomes \nBin size = ", window_size)) +
        xlab("Chromosome") + ylab("Log2 fold change") + theme(plot.title = element_text(size = 15)) + theme_cowplot()
      
      ggsave(paste0(log2fc_vector[log2fc_i], "_All_chr_violin.pdf"), width = 10, height = 7,  dpi = 500, path = paste0(BLISS_run, "/", window_size, "/Log2FC_analysis"))
      
    }
    
  }
  
}

