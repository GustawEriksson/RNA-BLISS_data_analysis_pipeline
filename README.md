README


                    ______ _   _   ___        ______ _     _____ _____ _____   
                    | ___ \ \ | | / _ \       | ___ \ |   |_   _/  ___/  ___|  
                    | |_/ /  \| |/ /_\ \______| |_/ / |     | | \ `--.\ `--.   
                    |    /| . ` ||  _  |______| ___ \ |     | |  `--. \`--. \  
                    | |\ \| |\  || | | |      | |_/ / |_____| |_/\__/ /\__/ /  
                    \_| \_\_| \_/\_| |_/      \____/\_____/\___/\____/\____/   


                    ______      _                            _           _     
                    |  _  \    | |                          | |         (_)    
                    | | | |__ _| |_ __ _    __ _ _ __   __ _| |_   _ ___ _ ___ 
                    | | | / _` | __/ _` |  / _` | '_ \ / _` | | | | / __| / __|
                    | |/ / (_| | || (_| | | (_| | | | | (_| | | |_| \__ \ \__ \
                    |___/ \__,_|\__\__,_|  \__,_|_| |_|\__,_|_|\__, |___/_|___/
                                                                __/ |          
                                                               |___/           
                    ______ _            _ _                                    
                    | ___ (_)          | (_)                                   
                    | |_/ /_ _ __   ___| |_ _ __   ___                         
                    |  __/| | '_ \ / _ \ | | '_ \ / _ \                        
                    | |   | | |_) |  __/ | | | | |  __/                        
                    \_|   |_| .__/ \___|_|_|_| |_|\___|                        
                            | |                                                
                            |_|                                                

                    A data analysis pipeline for RNA-sequencing and BLISS data

Author: Gustaw Eriksson

Date: 2020-05-14

Version: 1.0

Contact: gustaw.eriksson@ki.se

# Introduction:
The RNA-BLISS data analysis pipeline for RNA-sequencing and BLISS data, is used
to  process and analyse QoRTs formatted RNA-sequencing and BED-formatted BLISS
data. The pipeline processes the datasets separately, performing differential
gene expression analysis on the RNA-seq data and mapping chromosome-wide DNA 
double-strand breaks (DSB), before coupling the datasets. After coupling, the 
pipeline analysis and outputs information whether there is an relationship 
between gene expression and the accumulation of DSB.

In its current version, the pipeline requires a strict file directory
architecture which construction is mapped out below. To build it, requires the
use of Bedtools (v.2.29.2) and separate input in bash (v.3.2.57), also described
below.

The pipeline consists of several R-scripts, to be executed separately in an set
order to reproduce the computational pipeline described in "Understanding DNA
double-strand breaks and genome fragility across neurodifferentiation"
manuscript by Gustaw Eriksson (2020). For data to reproduce the computational 
pipeline, please read the Data availability section further down. The scripts where 
programmed using R (v.3.6.4) in RStudio (v.1.2.5033) with all packages installed 
with Bioconductor (v.3.10).   

The computational pipeline has only been tested on in-house produced data
derived from an in-vitro model of neurodifferentiation, involving three different 
groups (NES, progenitor, neuron) with three replicates per group. Running the pipeline 
on a different experimental layout involves risks and errors. Therefore, the current 
version is mere an experimental code NOT to be used without acknowledging the risks for 
errors and false results.

# Usage:
Bellow follows instruction on how to run the RNA-BLISS data analysis pipeline
and reproduce the results of "Understanding DNA double-strand breaks and genome
fragility across neurodifferentiation". It shows how the file directory was built, what 
the different folder contain and how the scripts were used:

  ## File directory architecture:
  The pipeline requires a set file directory architecture of three folders within
  the RNA-BLISS_Data_analysis_Pipeline folder:

### /Data/: 
The data folder should contain four folders called Bedtools_windows, BLISS_bed, 
hg19 and RNA_data:

  #### /Data/Bedtools_windows/: 
  Contains the Bedtools makewindow BED files further explained below. 
  The BED files should be of the different window sizes the user want to test and apply 
  on the BLISS data.

  The Bedtools (v.2.29.2) makewindow command was used in Bash in the following way, 
  where N is the window size:
    
      > bedtools makewindows -g hg19.genome -w N > N_window

  The generated windows for the manuscript were:
      
      BLISS = 10kb, 50kb, 100kb, 150kb, 200kb, 250kb, 300kb
      sBLISS = 1kb, 5kb, 10kb, 15kb, 20kb, 25kb, 30kb, 35kb, 40kb, 45kb, 50kb

  #### /Data/BLISS_bed/: 
  In the manuscript the program was written and applied on, the folder 
  contained two subfolders: /BLISS_bed/B138/ and /BLISS_bed/sBLISS/,containing the BLISS 
  BED files.

  #### /Data/hg19/: 
  Contained hg19.fa (full hg19 reference genome) and hg19.genome (chromosome sizes) 
  downloaded from USSC and gencode.v19.annotation.gtf (annotated hg19) downloaded from 
  GENCODE. The folder also contains a subfolder /annotation/ with the gencode.v19.annotation.gtf 
  converted to a .rds object called hg19_Gencode19_annotations.all.genes.regions.

  #### /Data/RNA_data/: 
  The folder is divided into two subfolders, RNA_BAM containing the RNA-sequencing 
  BAM files and RNA_qorts, containing the QoRTs formatted RNA-sequencing files.

### /Output/: 
When the program is first loaded, this folder should be empty. Running the R-scripts 
will generate subfolder within this folder containing output data and results.

### /Scripts/: 
The folder shall contain all the scripts part of the RNA-BLISS-Data_analysis_Pipeline. 
The script are divided between four different folders:

  #### /Scripts/BLISS_analysis/: 
  Contains the Automate_BLISS_QC.R and BLISS_functions.R. The 
  Automate_BLISS_QC.R processes and run quality control on the BLISS data. The code is to be run 
  on both the BLISS and sBLISS data, which is default. The code calls functions from the 
  BLISS_functions.R script. When Automate_BLISS_QC.R is open in RStudio, the user can change 
  input and output directories, and which window sizes it is going to use as input when processing 
  the BLISS data. It is required that these window sizes have been generated earlier in the 
  Bedtools makewindows step.

  Depending on the result from the Pearson and Spearman correlation, one window size is chosen for 
  downstream analysis. The default window for BLISS and sBLISS in the current version are 150kb 
  and 10kb respectively.

  #### /Scripts/RNA_analysis/: 
  Within the folder, are two scripts, the Automate_RNA_analysis.R which runs the 
  DESeq2 (v1.24.0) pipeline including quality control and GSEA by clusterProfiler (v.3.0.4), and 
  RNA_analysis_tools.R which keeps functions used in the Automate_RNA_analysis.R script.

  When Automate_RNA_analysis.R is used, the default thresholds are:
    
    log2fc = -1/1
    p-value = <0.05
    TPM = 1

  The default data transformation method is variance-stabilizing transformation (VST). When 
  running the pipeline on a different set of data, it is highly recommended to inspect the data 
  transformation QC to determine whether one of the two other data transformation methods are 
  better: regularised log transformation (rlog) or normalised count transformation 
  (normTransform), besides VST.

  After setting data transformation method, study the QC of the samples, which are the Euclidean 
  distances and PCA plot, so that the samples are of suitable quality.

  The code generates figures, plots and tables illustrating differential gene expression of genes. 
  Genes are select after the gene expression pattern and if its differential exoression is above 
  or below the log2fc threshold and has an significant p-value, across the different cell stages. 
  One of the generated heatmaps consists of pre-selected genes that have been reported to be 
  fragile and accumulate DSBs by Wang et al. (2020). This pre-selected gene list can be changed by 
  the user to show the differential gene expression of other genes.

  #### /Scripts/DSB_mapping/: 
  Both Automate_BLISS_QC.R and Automate_RNA_analysis.R has to have been executed 
  before running the scripts found in /DSB_mapping/. chromosome_wide_DSB.R generates output 
  showing the DSB frequency across chromosome 1-22 and X, whilst also calculating the log2 fold 
  change of DSB frequency between bins of set window size, between the different experimental 
  groups. The DSB_mapping_tools.R keeps function used in chromosome_wide_DSB.R. The script also 
  depend on /Data_handling/core_functions.R which is described below.

  Function_DSB_TSS.R merges data from the RNA-sequencing and BLISS data analysis to map the 
  accumulation of DSB in and surrounding the transcription start site (TSS) of genes and is 
  dependent of Subread package (v.2.0.1). The code also analyses the overlap of extracted genes 
  from the RNA-sequencing data analys part of the script with genes found within bins having 
  increased or decreased frequency of DSB between cell stages. For a bin to be determined has 
  having an increased or decreased frequency of DSB, it must have a log2fc below or above the set 
  log2fc threshold, which default setting in Function_DSB_TSS.R is -1.5/1.5.

  #### /Scripts/Data_handling/: 
  The folder contains the core_functions.R by Jesko Wagner, that loads BLISS-data 
  and is used in /DSB_mapping_tools.R.

## Data availability:
The data used to program this data analysis pipeline, was generated at the Department of Medical 
Biochemistry and Biophysics, Karolinska Institute by the Laboratory of Quantative Genome Biology led by 
Nicola Crosetto. The data is currently not public.

The author of the program is Gustaw Eriksson. You can reach me at gustaw.eriksson@ki.se if you have any 
questions regarding the pipeline or anything else.
