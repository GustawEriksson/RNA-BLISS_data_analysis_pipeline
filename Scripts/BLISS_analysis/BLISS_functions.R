# Author: Gustaw Eriksson
# Date: 2020-05-14

# Description: Function used in Automate_BLISS_QC.R script doing qualtiy control of 
# BLISS and sBLISS libraries to eachother.

library(data.table)
library(stringr)
library(plyr)
library(dplyr)
library(tidyr)
library(tibble)
library("Hmisc")
library("GGally")
library(GenomicRanges)
library(bedr)
library(HelloRanges)
library(naturalsort)
library(cowplot)

statistics_sample_p_value <- function(input_DF) {
  data_rcorr = as.matrix(input_DF)
  input_DF_rcorr = rcorr(data_rcorr)
  input_sample_p_value = input_DF_rcorr[["P"]]
  
  return(input_sample_p_value)
}

statistics_spearman <- function(input_DF, window_size) {
  spearman_ggcorr = ggcorr(input_DF, method = c("pairwise", "spearman"), midpoint = NULL, label = TRUE, label_round = 2, hjust = 1, size = 4, layout.exp = 4)
  spearman_ggcorr + ggtitle(paste0("Spearman correlation ", window_size, " windows"))
  return(spearman_ggcorr)
}

statistics_pearson <- function(input_DF, window_size) {
  pearson_ggcorr = ggcorr(input_DF, method = c("pairwise", "pearson"), midpoint = NULL, label = TRUE, label_round = 2, hjust = 1, size = 4, layout.exp = 4)
  pearson_ggcorr + ggtitle(paste0("Pearson correlation ", window_size, " windows"))
  return(pearson_ggcorr)
}

coefficient_of_variation <- function(input_DF, add_to_DF = TRUE) {
  co.var <- function(x) (100*sd(x)/mean(x))
  input_DF_CV = apply(input_DF, 1, co.var)
  return(input_DF_CV)
  
}

plot_Correlation_ggpair = function(x, suffix, w = 15, h = 15) {
  continuousCorrelation_wLine = function(data, mapping) {
    return(ggplot(data, mapping
    ) + geom_point(size = .5, alpha = 0.4
    ) + geom_abline(intercept = 0, slope = 1,
                    color = "red", linetype = "dashed"))}
  
  continuousCorrelation_mixed = function(data, mapping) {
    Xcol = rlang::get_expr(mapping$x)
    x = data[[Xcol]]
    
    Ycol = rlang::get_expr(mapping$y)
    y = data[[Ycol]]
    
    pcor = cor(x, y, use = "complete.obs", method = "pearson")
    scor = cor(x, y, use = "complete.obs", method = "spearman")
    
    return(ggplot(data, mapping) + geom_blank(
    ) + geom_text(
      x = mean(c(min(x, na.rm = T), max(x, na.rm = T))),
      y = mean(c(min(y, na.rm = T), max(y, na.rm = T))),
      label = sprintf("Pearson: %.3f\nSpearman: %.3f", pcor, scor),
      family = "Helvetica", size = 5, color = "black",
      fontface = "plain"))
  }
  
  p = ggpairs(x, columns = 1:ncol(x), title = paste0
              ("Pearson and Spearman correlation for ", suffix),
              lower = list(continuous = continuousCorrelation_wLine),
              upper = list(continuous = continuousCorrelation_mixed),
              axisLabels = "show", columnLabels = names(x)[1:ncol(x)]
  ) + theme_minimal()
}

generate_corr_matrix <- function(input, suffix) {
  corr_matrix = as.data.table(input, keep.rownames = "window") %>% rename("window" = "rn") %>% column_to_rownames(var = "window") %>% 
    dplyr::select(matches(suffix)) %>% rownames_to_column(var = "window") %>% filter(grepl(suffix, window))
  return(corr_matrix)
}

plot_CV <- function(CV, window_size, suffix) {
  
  n_NA = paste0("Number of NA = ", sum(is.na(CV)))
  CV_df = as.data.frame(CV) %>% rownames_to_column(., var = "window") %>% setnames(., c("window", "CV")) %>% na.omit(.) %>% dplyr::arrange(desc(CV)) %>%
    mutate(window=factor(window, levels = window))
  CV_plot = ggplot(CV_df, aes(x = window, y = CV)) + geom_point(colour = "red") + ylim(0, 400) + ylab("CV%") + xlab("Ordered bins by descending CV%") +
    ggtitle(label = paste0("CV% for ", window_size, " bins in ", suffix), subtitle = n_NA)
  
}

plot_Correlation_matrix <- function(indir, window_vector, input_matrix, 
                                    sBLISS = FALSE, BLISS_suffix) {
  
  length_window_vector = length(window_vector)
  n_windows = 0
  matrix_value_list = list()
  load_matrix_list = list()
  
  sBLISS_100_8 = c()
  sBLISS_100_12 = c()
  sBLISS_300_8 = c()
  sBLISS_300_12 = c()
  
  BLISS_1_2 = c()
  BLISS_1_3 = c()
  BLISS_2_3 = c()
  
  # Load in matrix in loop over window vector
  for (window in window_vector) {
    n_windows = n_windows + 1
    matrix_input = paste0(indir, window, "/")
    matrix_target = list.files(matrix_input, pattern = paste0(input_matrix))
    matrix_input = paste0(matrix_input, matrix_target)
    target_matrix = fread(file = matrix_input)
    window_flag = names(target_matrix[,1]) == "window"
    if (window_flag == TRUE) {
      target_matrix = target_matrix %>% column_to_rownames(., var = "window") 
    }
    
    if (sBLISS == TRUE) {
      sBLISS_100_8 = c(sBLISS_100_8, target_matrix[2,1])
      sBLISS_100_12 = c(sBLISS_100_12, target_matrix[4,3])
      sBLISS_300_8 = c(sBLISS_300_8, target_matrix[6,5])
      sBLISS_300_12 = c(sBLISS_300_12, target_matrix[8,7])
      
    } else if (sBLISS == FALSE) {
      BLISS_1_2 = c(BLISS_1_2, target_matrix[2,1])
      BLISS_1_3 = c(BLISS_1_3, target_matrix[3,1])
      BLISS_2_3 = c(BLISS_2_3, target_matrix[2,3])
    }
  }
  
  if (sBLISS == TRUE) {
    plot_df = data.frame(BLISS=rep("sBLISS", 40),replicate=rep(c("sBLISS_100_8", "sBLISS_100_12", "sBLISS_300_8", "sBLISS_300_12"), each=length_window_vector), 
                         window_size = rep(window_vector, 4), correlation = c(sBLISS_100_8, sBLISS_100_12, sBLISS_300_8, sBLISS_300_12))
  } else if (sBLISS == FALSE) {
    plot_df = data.frame(BLISS=rep(BLISS_suffix, 30), replicate=rep(c(paste0(BLISS_suffix,"_1_vs._2"), paste0(BLISS_suffix,"_1_vs._3"), paste(BLISS_suffix,"_2_vs._3")), each=length_window_vector), 
                         window_size = rep(window_vector, 3), correlation = c(BLISS_1_2, BLISS_1_3, BLISS_2_3))
  }
  
  return(plot_df)
  #ggplot(plot_df, aes(x = factor(window_size, level = window_vector), y=correlation, group=replicate)) + 
  #  geom_line(aes(color=replicate)) + geom_point(aes(color=replicate)) + theme_cowplot()
  
}

bedtools_window_merged_count <- function(bed_indir, bed_window_dir, bed_pattern = ".bed$", exclude_chrY = TRUE, 
                                         order_chr = TRUE, sep_window_col = FALSE) {
  
  ## Author: Gustaw Eriksson
  ## Date: 2019-11-11
  ## Description: Input directory of bed files with DSB count column and path to bin window bed.
  ## The function pairs the bed files to the bin window file and merges the files to one matrix 
  ## containing window and number of DSB in the window.
  
  bed_window = fread(bed_window_dir)
  
  if (order_chr == TRUE){
    bed_window = bed_window[naturalorder((bed_window$V1), decreasing = FALSE),]
  }
  
  bed_indir_list = list.files(bed_indir, pattern = bed_pattern)
  bed_files_list = list()
  
  for (bed in 1:length(bed_indir_list)) {
    bed_key = bed_indir_list[[bed]]
    bed_file = fread(file = paste0(bed_indir, "/", bed_indir_list[bed]))
    bed_files_list[[bed_key]] = bed_file
    
  }
  
  colnames(bed_window) = c("chr", "start", "end")
  bed_window_GR = makeGRangesFromDataFrame(bed_window)
  
  first_bed = TRUE

  for (bed in 1:length(bed_files_list)) {
    bed_df = bed_files_list[[bed]]
    colnames(bed_df) = c("V1", "V2", "V3", "V4")
    bed_df = mutate_at(bed_df, vars(V1), as.character) %>% filter(., grepl('chr', V1)) %>% filter(., !grepl('chrM', V1))
    
    if (exclude_chrY == TRUE) {
      bed_df = bed_df %>% filter(., !grepl('chrY', V1))
    }
    
    if (first_bed == TRUE) {
      first_bed = FALSE
      
      main_overlap_df = as.data.table(bed_window_GR) %>% set(., , 4:5, NULL) %>% unite("window", start, end, sep='-', remove = TRUE) %>%
        unite("window", seqnames, window, sep= ' ', remove = TRUE)
      
    } 
    
    bed_GR = makeGRangesFromDataFrame(bed_df, keep.extra.columns = T, seqnames.field = "V1", start.field = "V2", end.field = "V3")
    overlap_GR <-findOverlapPairs(bed_window_GR, bed_GR)
    
    append_overlap_df = as.data.table(overlap_GR) %>% set(., , 4:11, NULL) %>% unite("window", first.start, first.end, sep='-', remove = TRUE) %>%
      unite("window", first.seqnames, window, sep= ' ', remove = TRUE) %>% group_by(window, .drop = FALSE) %>% dplyr::summarise(second.V4 = sum(second.V4))
    
    bed_name = names(bed_files_list[bed]) %>% str_extract(., "[^.]+")
    colnames(append_overlap_df) = c("window", bed_name)
    main_overlap_df = left_join(main_overlap_df, append_overlap_df, by = "window", copy = FALSE) %>% replace(is.na(.), 0)
    
  }
  
  if (exclude_chrY == TRUE) {
    main_overlap_df = main_overlap_df %>% filter(., !grepl('chrY', window))
  }
  
  if (sep_window_col == FALSE) {
    return(main_overlap_df)
  } else if (sep_window_col == TRUE) {
    main_overlap_df = separate(main_overlap_df, window, c("Chr", "Start", "End"), sep = "[ ]|[-]")
    main_overlap_df$Chr = gsub('chr', '', main_overlap_df$Chr)
    return(main_overlap_df)
  }
  
}

filter_CV_rank_bin <- function(normalised_bin_DT, target_BLISS_pattern, CV_threshold) {
  
  library('data.table')
  library('dplyr')
  library('tibble')
  
  # Subset selected BLISS run
  normalised_bin_DT = setDT(normalised_bin_DT, keep.rownames = TRUE)
  subset_DT = normalised_bin_DT %>% column_to_rownames(., var = "window") %>% dplyr::select(matches(target_BLISS_pattern))

  # Running CV% statistics on data table
  subset_CV = coefficient_of_variation(subset_DT)

  # Filtering CV% vector with set CV% threshold
  subset_CV = as.data.table(subset_CV, keep.rownames = TRUE) %>% rename(., "CV" = "subset_CV") %>% rename(., "window" = "rn") %>% 
    filter(., CV <= CV_threshold)
  
  # Match filtered CV% vector against normalised BLISS DT to get normalised counts with CV% within threshold
  mean_subset_DT = subset_DT %>% rownames_to_column(., var = "window") %>% left_join(., subset_CV, by = "window") %>%
    filter(., !is.na(CV)) %>% dplyr::select(-CV) %>% setDT(.) %>% .[, .(Mean = rowMeans(.SD)), by = window]
  
  # Run mean DSB count and l join in normalised BLISS DT
  normalised_filtered_bin_DT = subset_DT %>% rownames_to_column(., var = "window") %>% 
    inner_join(., mean_subset_DT, by = "window") %>% setDT(.)
  
  return(normalised_filtered_bin_DT)
  
}

# When doing this control on the sBLISS, first compare the libraries to each other as they
# have been produced using different volumes of chemicals and number of PCR cycles. 
LOAD_NORMALISED_BIN <- function(indir, outdir, window_size, new_window = FALSE, subset_DF = NULL, bed_pattern = ".bed$") {
  
  window_file = paste0("/RNA-BLISS_Data_analysis_Pipeline/Data/Bedtools_windows/", window_size,"_window")

  if (new_window == TRUE) {
    ## Send list of bed files to bedtools function
    bin_DT = bedtools_window_merged_count(bed_indir = indir, bed_window_dir = window_file, bed_pattern = bed_pattern)
    
    ## Diving each DSB count with total DSB count
    normalised_bin_DT = bin_DT %>% mutate_if(., is.numeric, funs(((./sum(.))*1000000)))
    
    ## Saving the data.table as csv
    print(paste0(outdir, "   OUTDIR FOR BIN_DT"))
    fwrite(x = bin_DT, file = paste0(outdir, "/", bed_pattern,"_", window_size,"_bins_raw_count.csv"), sep = '\t')
    fwrite(x = normalised_bin_DT, file = paste0(outdir, "/", bed_pattern,"_", window_size,"_bins_normalised_count.csv"), sep = '\t')
    
    
  } else if (new_window == FALSE) {
    raw_count_input = list.files(indir, pattern = paste0(window_size,"_bins_raw_count.csv"))
    bin_DT = fread(file = paste0(indir, raw_count_input))
    bin_DT = bin_DT[order(bin_DT$window),]
    normalised_count_input = list.files(indir, pattern = paste0(window_size,"_bins_normalised_count.csv"))
    normalised_bin_DT = fread(file = paste0(indir,normalised_count_input))
    
  }
  
  ## Setting window column as rowname
  bin_DT = bin_DT %>% column_to_rownames(., var = "window")
  normalised_bin_DT = normalised_bin_DT %>% column_to_rownames(., var = "window")
  
  ## Creating heatmap with only specific libraries

  if (is.null(subset_DF) == FALSE) {
    
    normalised_bin_DT = subset(normalised_bin_DT, select = -c(subset_DF))
    
  }
  
  return(normalised_bin_DT)
  
}

BLISS_BIN_QC <- function(normalised_bin, outdir, window_size, suffix, new_window = FALSE, CV_plot = TRUE, gg_pairs_plot = TRUE, 
                         CV_0_plot = TRUE, gg_pairs_0_plot = TRUE) {
  
  # Running sample_p_value on DT. Are samples random and independent
  target_normalised_p_value = statistics_sample_p_value(normalised_bin)
  fwrite(target_normalised_p_value, paste0(outdir, "/", window_size, "_sample_p_value.csv"), sep = '\t')
  
  # Spearman correlation block
  spearman_matrix = cor(normalised_bin, method = "spearman")
  fwrite(spearman_matrix, file = paste0(outdir, "/", window_size,"_spearman_matrix.csv"), sep = '\t')
  spearman_target_normalised_ggcorr = statistics_spearman(normalised_bin, window_size)
  ggsave(file = paste0(outdir, "/", window_size, '_spearman_correlation.png'), plot = spearman_target_normalised_ggcorr, width=8, height=10, units="in")

  # Pearson correlation block
  pearson_matrix = cor(normalised_bin, method = "pearson")
  fwrite(pearson_matrix, file = paste0(outdir, "/", window_size,"_pearson_matrix.csv"), sep = '\t')
  pearson_target_normalised_ggcorr = statistics_pearson(normalised_bin, window_size)
  ggsave(file = paste0(outdir, "/", window_size, '_pearson_correlation.png'), plot = pearson_target_normalised_ggcorr, width=8, height=10, units="in")
  
  target_DF_no_0_list = list()
  
  for (target in suffix) {
    
    if (length(suffix) > 1) {
      # Subsetting normalised bin DT
      target_DF = normalised_bin %>% dplyr::select(contains(target))
      
      # Line plots for pearson and spearman correlation
      target_pearson_matrix = generate_corr_matrix(input = pearson_matrix, suffix = target)
      fwrite(target_pearson_matrix, paste0(outdir, "/", window_size,"_", target, "_pearson_matrix.csv"))
      target_spearman_matrix = generate_corr_matrix(input = spearman_matrix, suffix = target)
      fwrite(target_spearman_matrix, paste0(outdir, "/", window_size,"_", target, "_spearman_matrix.csv"))
      
    } else if (length(suffix) == 1) {
      target_DF = normalised_bin
    }
    
    ##Go through each row and determine if a value is zero
    row_sub = apply(target_DF, 1, function(row) all(row !=0 ))
    ##Subset as usual
    target_DF_no_0 = target_DF[row_sub,]
    
    # Pearson correlation block
    target_DF_no_0_pearson_matrix = cor(target_DF_no_0, method = "pearson")
    fwrite(target_DF_no_0_pearson_matrix, file = paste0(outdir, "/", window_size,"_", target, "_pearson_no_0_matrix.csv"))
    target_DF_no_0_spearman_matrix = cor(target_DF_no_0, method = "spearman")
    fwrite(target_DF_no_0_spearman_matrix, file = paste0(outdir, "/", window_size,"_", target, "_spearman_no_0_matrix.csv"))
    
    if (gg_pairs_plot == TRUE){
      #Running ggpairs
      target_ggpairs = plot_Correlation_ggpair(target_DF, suffix = target)
      ggsave(file = paste0(outdir, "/", window_size, '_', target, '_ggpairs.png'), plot = target_ggpairs, width=22, height=22, units="in")
      
    }
    
    if (CV_plot == TRUE){
      
      target_CV= coefficient_of_variation(target_DF)
      target_CV_plot = plot_CV(CV = target_CV, window_size = window_size, suffix = target)
      ggsave(filename = paste0(outdir, "/", window_size, '_plot_percantage_CV_', target,'.png'), plot = target_CV_plot, width=8, height=10, units="in")
      
    }
    
    if (CV_0_plot == TRUE | gg_pairs_0_plot == TRUE) {
      ##Go through each row and determine if a value is zero
      row_sub = apply(target_DF, 1, function(row) all(row !=0 ))
      ##Subset as usual
      target_DF_no_0 = target_DF[row_sub,]
      
      # Running correlation with no zero
      spearman_target_DF_no_0 = statistics_spearman(target_DF_no_0, window_size)
      ggsave(file = paste0(outdir, "/", window_size, '_', target, '_spearman_correlation_no_0.png'), plot = spearman_target_DF_no_0, width=8, height=10, units="in")
      pearson_target_DF_no_0 = statistics_pearson(target_DF_no_0, window_size)
      ggsave(file = paste0(outdir, "/", window_size, '_', target, '_pearson_correlation_no_0.png'), plot = pearson_target_DF_no_0, width=8, height=10, units="in")
      
      
      if (gg_pairs_0_plot == TRUE){
        target_ggpairs = plot_Correlation_ggpair(target_DF_no_0, suffix = target)
        ggsave(file = paste0(outdir, "/", window_size, '_', target, '_ggpairs_no_0.png'), plot = target_ggpairs, width=22, height=22, units="in")
        
      }
      
      if (CV_0_plot == TRUE) {
        target_CV = coefficient_of_variation(target_DF_no_0)
        target_CV_plot = plot_CV(CV = target_CV, window_size = window_size, suffix = target)
        ggsave(filename = paste0(outdir, "/", window_size, '_', target,'_plot_percantage_CV_no_0.png'), plot = target_CV_plot, width=8, height=10, units="in")
      }
    }
  }
}

EXTRACT_TARGET_MATRIX <- function(indir, input_matrix) {

  matrix_target = list.files(indir, pattern = paste0(input_matrix))
  matrix_input = paste0(indir, matrix_target)
  target_matrix = fread(file = matrix_input)
  if ((names(target_matrix[,1]) == "window") == TRUE) {
    target_matrix = target_matrix %>% column_to_rownames(., var = "window") 
  }
  return(target_matrix)
}

BLISS_window_correlation_DF = function(indir, windows, Cell_type, input_matrix) {
  first_matrix = TRUE
  n_x = 1
  for (window in windows) {
    print(window)
    indir_window = paste0(indir, window, "/", Cell_type,"_QC/")
    print(indir_window)
    window_matrix = as.matrix(EXTRACT_TARGET_MATRIX(indir = indir_window, input_matrix = input_matrix))
    if (first_matrix == TRUE) {
      first_matrix = FALSE
      ncol_matrix = ncol(window_matrix)
      len_windows = length(windows)
      if (ncol_matrix == 3) {
        matrix_pos = c("1-2", "1-3", "2-3")
        x_corr = c(rep(1, 2), 2)
        y_corr = c(2, rep(3, 2))
        corr_DF = data.frame(Cell_type = rep(Cell_type, (ncol_matrix*len_windows)), Replicate = rep(matrix_pos, len_windows), 
                             window_size = rep(windows,each = ncol_matrix), correlation =rep(NA, len_windows))
      } else if (ncol_matrix == 4) {
        matrix_pos = c("1-2", "1-3", "1-4", "2-3", "2-4", "3-4")
        x_corr = c(rep(1, 3), rep(2, 2), 3)
        y_corr = c(2, 3, 4, 3, 4, 4)
        corr_DF = data.frame(Cell_type = rep(Cell_type, ((ncol_matrix*1.5)*len_windows)), Replicate = rep(matrix_pos, len_windows), 
                             window_size = rep(windows,each = ncol_matrix*1.5), correlation =rep(NA, len_windows))
      } else if (ncol_matrix == 5) {
        matrix_pos = c("1-2", "1-3", "1-4", "1-5", "2-3", "2-4", "2-5", "3-4", "3-5", "4-5")
        x_corr = c(rep(1, 4), rep(2, 3), rep(3, 2), 4)
        y_corr = c(2, 3, 4, 5, 3, 4, 5, 4, 5, 5)
        corr_DF = data.frame(Cell_type = rep(Cell_type, ((ncol_matrix*2)*len_windows)), Replicate = rep(matrix_pos, len_windows), 
                             window_size = rep(windows,each = ncol_matrix*2), correlation =rep(NA, len_windows))
      } else if (ncol_matrix == 8) {
        print("MATRIX RUNNING")
        matrix_pos = c("1-2", "1-3", "1-4", "1-5", "1-6", "1-7", "1-8", 
                       "2-3", "2-4", "2-5", "2-6", "2-7", "2-8", 
                       "3-4", "3-5", "3-6", "3-7", "3-8", 
                       "4-5", "4-6", "4-7", "4-8",
                       "5-6", "5-7", "5-8", 
                       "6-7", "6-8", 
                       "7-8")
        x_corr = c(rep(1, 7), rep(2, 6), rep(3, 5), rep(4,4), rep(5,3), rep(6, 2), 7)
        y_corr = c(2, 3, 4, 5, 6, 7, 8, 3, 4, 5, 6, 7, 8, 4, 5, 6, 7, 8, 5, 6, 7, 8, 6, 7, 8, 7, 8, 8)
        corr_DF = data.frame(Cell_type = rep(Cell_type, ((ncol_matrix*3.5)*len_windows)), Replicate = rep(matrix_pos, len_windows), 
                             window_size = rep(windows,each = ncol_matrix*3.5), correlation =rep(NA, len_windows))
      }
    }
    for (i in 1:length(x_corr)) {
      corr_value = window_matrix[x_corr[i], y_corr[i]]
      corr_DF[n_x, 4] = corr_value
      n_x = n_x + 1
    }
  }
  return(corr_DF)
}

PLOTTING_WINDOW_CORR <- function(indir = "/RNA-BLISS_Data_analysis_Pipeline/Output/BLISS/", 
                                 Cell_type = c("NES", "Progenitor", "Neuron"), BLISS_run = "B138", 
                                 windows = c("500Kb", "600Kb", "700Kb", "800Kb", "900Kb", "1Mb"), 
                                 input_matrix = "pearson_matrix.csv", plotting = TRUE) {
  
  if (grepl("pearson", input_matrix) == TRUE) {
    print("WORKING")
    corr_type = "Pearson"
  } else if (grepl("spearman", input_matrix) == TRUE) {
    corr_type = "Spearman"
  }
  
  first_cell = TRUE
  for (cell in Cell_type) {
    matrix_indir = paste0(indir, BLISS_run, "/")
    print(matrix_indir)
    cell_matrix = BLISS_window_correlation_DF(indir = matrix_indir, windows = windows, Cell_type = cell, input_matrix = input_matrix)
    
    if (plotting == TRUE){
      cell_matrix$window_size = gsub("[^0-9]", "", cell_matrix$window_size)
      window_levels = cell_matrix$window_size[!duplicated(cell_matrix$window_size)]
      corr_plot_cell = ggplot(cell_matrix, aes(x = factor(window_size, level = window_levels), y=correlation, group=interaction(Replicate, Cell_type))) + 
        geom_line(aes(color=Replicate, linetype=Cell_type)) + 
        geom_point(aes(color=Replicate, shape=Cell_type)) + theme_cowplot() +
        labs(title=paste0(corr_type," correlation of different bin sizes between \nreplicates in ", cell),
             x="Bin window size (Kb)", y = paste0(corr_type," correlation")) + theme(axis.text.x=element_text(angle=45, vjust = 0.80))
      ggsave(filename = paste0(matrix_indir, corr_type, "_", cell, "_correlation_plot.pdf"), plot = corr_plot_cell, 
             device = "pdf")
    }
    
    
    if (first_cell == TRUE) {
      first_cell = FALSE
      plot_df = cell_matrix
    } else if (first_cell == FALSE) {
      plot_df = rbind(plot_df, cell_matrix)
    }
    
  }
  
  if (plotting == TRUE){
    plot_df$window_size = gsub("[^0-9]", "", plot_df$window_size)
    corr_plot_cell = ggplot(plot_df, aes(x = factor(window_size, level = window_levels), y=correlation, group=interaction(Replicate, Cell_type))) + 
      geom_line(aes(color=Replicate, linetype=Cell_type)) + geom_point(aes(color=Replicate, shape=Cell_type)) + theme_cowplot() +
      labs(title=paste0(corr_type," correlation of different bin sizes between \nreplicates across neurodifferentiation"),
           x="Bin window size (Kb)", y = paste0(corr_type," correlation")) + theme(axis.text.x=element_text(angle=45, vjust = 0.8)) +
      geom_hline(yintercept = 0.8, linetype = "dashed")
    ggsave(filename = paste0(indir, corr_type, "_", BLISS_run, "_all_cell_types.pdf"), plot = corr_plot_cell, device = "pdf")
  }
  
  return(plot_df)
}

