
BLISS_windows_filt = c("10Kb", "50Kb", "100Kb", "150Kb", "200Kb", "250Kb", "300Kb")

sBLISSwindows_filt = c("1Kb", "5Kb", "10Kb", "15Kb", "20Kb", "25Kb", "30Kb", "35Kb", "40Kb", "45Kb", "50Kb")

sBLISS_run = run_QC_BLISS(BLISS_run = "sBLISS", cell_type_vector = c("NES"), 
                          #windows = c("10Kb", "20Kb", "30Kb", "40Kb", "50Kb", "60Kb", "70Kb", "80Kb", "90Kb", "100Kb", "200Kb", "300Kb", "400Kb", "500Kb"),
                          windows = sBLISS_windows_filt,
                          indir = "/RNA-BLISS_Data_analysis_Pipeline/Data/BLISS_bed/", 
                          outdir = "/RNA-BLISS_Data_analysis_Pipeline/Output/")

pearson_sBLISS_plots = PLOTTING_WINDOW_CORR(indir = "/RNA-BLISS_Data_analysis_Pipeline/Output/BLISS/",
                                            Cell_type = "NES", BLISS_run = "sBLISS",
                                            windows = sBLISS_windows_filt, input_matrix = "pearson_matrix.csv")
spearman_sBLISS_plots = PLOTTING_WINDOW_CORR(indir = "/RNA-BLISS_Data_analysis_Pipeline/Output/BLISS/",
                                             Cell_type = "NES", BLISS_run = "sBLISS",
                                             windows = sBLISS_windows_filt, input_matrix = "spearman_matrix.csv")

BLISS_RUN = run_QC_BLISS(windows = BLISS_windows_filt, indir = "/RNA-BLISS_Data_analysis_Pipeline/Data/BLISS_bed/", 
                         outdir = "/RNA-BLISS_Data_analysis_Pipeline/Output/")

pearson_BLISS_plots = PLOTTING_WINDOW_CORR(windows = BLISS_windows_filt, input_matrix = "pearson_matrix.csv")
spearman_BLISS_plots = PLOTTING_WINDOW_CORR(windows = BLISS_windows_filt, input_matrix = "spearman_matrix.csv")


run_QC_BLISS <- function(BLISS_run, cell_type_vector = c("NES", "Progenitor", "Neuron"), 
                         windows, 
                         indir = "/RNA-BLISS_Data_analysis_Pipeline/Data/BLISS_bed/", 
                         outdir = "/RNA-BLISS_Data_analysis_Pipeline/Output/") {
  
  #source("/Users/gustaweriksson/Dokument/Bienko-Crosetto/Final_pipeline/Scripts/Quality_control/BLISS/BLISS_BIN_QC.2.R")
  source("/RNA-BLISS_Data_analysis_Pipeline/Scripts/BLISS_analysis/BLISS_functions.R")
  
  indir = paste0(indir, BLISS_run)
  
  setwd(outdir)
  dir.create("BLISS", showWarnings = FALSE)
  dir.create(paste0("BLISS/", BLISS_run), showWarnings = FALSE)
  #dir.create(paste0("BLISS/", BLISS_run, "/", window), showWarnings = FALSE)
  
  for (type in cell_type_vector) {
    print(type)
    for (window in windows) {
      
      outdir_tables = paste0("BLISS/", BLISS_run, "/", window, "/Count_tables")
      outdir_QC = paste0(outdir, "BLISS/", BLISS_run, "/", window, "/", type, "_QC")
      
      dir.create(paste0("BLISS/", BLISS_run, "/", window), showWarnings = FALSE)
      dir.create(outdir_tables, showWarnings = FALSE)
      dir.create(outdir_QC, showWarnings = FALSE)
      
      print(window)
      print(indir)
      if (grepl("sBLISS", BLISS_run) == FALSE) {
        normalised_bin = LOAD_NORMALISED_BIN(indir = indir, outdir = outdir_tables, 
                                             window_size = window, new_window = TRUE, subset_DF = NULL, bed_pattern = type)
      } else if (grepl("sBLISS", BLISS_run) == TRUE) {
        normalised_bin = LOAD_NORMALISED_BIN(indir = indir, outdir = outdir_tables, #outdir = outdir_window,
                                             window_size = window, new_window = TRUE, subset_DF = NULL)
      }
      
      print("Normalised bin done")
      run_QC = BLISS_BIN_QC(normalised_bin = normalised_bin, outdir = outdir_QC, #outdir = outdir_window,
                            window_size = window, suffix = BLISS_run, CV_plot = FALSE, gg_pairs_0_plot = FALSE, CV_0_plot = FALSE, gg_pairs_plot = FALSE)
      print("QC done")
    }
    
  }
  
  
}
