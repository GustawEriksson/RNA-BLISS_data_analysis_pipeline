### Project: Neuronal BLISS
#
### Aim: Have robust functions for the reading and processing of bed files for high-throughtput and reproducbility.
#
### Description: These functions allow for easy and reproducible
### loading and processing of data contained in bed files.
#
### Disclaimer: I use "groups" for individual experiments, "sample" for different sampels (e.g. NES, Progenitor, ...)
### and "replicate" or "rep" for distinguishing replicates of the same sampe type.
# 
### Example usage:
## Example 1: one group
# indir = "~/Desktop/user_folders/Jesko/bed/BLISS/B187/" # path containing folders for groups
# sep = "_"                             # what are sample description and replicate number separated by?
# annotation_file = "~/Desktop/user_folders/Jesko/annotations/hg19/general/hg19_Gencode19_annotations.all.genes.regions.rds" # path to GRanges of annotations
# 
# # load in all bed files into a nested, named list
# dt_list = load_bed_from_dir(indir, annotation_file = annotation_file, sep = sep)
# 
# # convert the nested list into a data.table that can be used for processing
# B187_dt = make_DT_from_list(dt_list)
#
## Example 2: multiple groups
# indir = ~/Desktop/user_folders/Jesko/bed/BLISS/ 
# multiple_groups = T                   # load in multiple groups
# annotation_file = "~/Desktop/user_folders/Jesko/annotations/hg19/general/hg19_Gencode19_annotations.all.genes.regions.rds" # path to GRanges of annotations
# 
# dt_list = load_bed_from_groups(indir = indir, multiple_groups = T, annotation_file = annotation_file)
# grouped_dt = make_DT_from_list(dt_list)

## Example 3: DT usage
## If wanted, the resulting data table can be used for plotting, e.g. pie charts
# library("ggplot2")
# p1 = grouped_dt %>%
#   .[,.(count = sum(libgenadj)), by=.(sample,annotation)] %>%
#   .[,.(annotation = annotation, percent=count/sum(count)*100), by=.(sample)] %>%
#   ggplot(aes(x="", y=percent, fill=annotation)) +
#   geom_bar(stat="identity") +
#   coord_polar("y", start=0) +
#   facet_wrap( ~ sample) +
#   geom_text(aes(label = paste0(round(percent,2), "%")), position = position_stack(vjust = 0.5))


### Depdencies
suppressMessages(library("tools"))
suppressMessages(library("rtracklayer"))
suppressMessages(library("purrr"))
suppressMessages(library("data.table"))
suppressMessages(library("GenomicRanges"))
suppressMessages(library("dplyr"))         # used exclusively for piping, so magrittr may be used instead
suppressMessages(library("parallel"))


### Aim: list_depth is a function to loop through a list and determine how many times it is nested.
### In other functions, it is used to determine how many groups (one or more) we are looking at.
list_depth = function(list) {
  if(! "list" %in% class(list)) return(0L)
  if(!length(list)) return(0L) # if the list is empty do not count it as nested
  return(1L+list_depth(list[[1]]))
}


### Aim: loading in bed files from a directory (indir), where the files have a group name and a replicate
### number (separated by (sep)). A group order can be given as character vector. Pre-selection of 
### chromosomes is possible, and, if libraryd, an integer of threads (p) to be used for fread.
### Example usage: load_bed_from_dir(indir) [where indir contains files of pattern "NES_1.bed"]

load_bed_from_dir = function(indir, sep = "_", samples = NULL, verbose = F, p = NULL) {
  
  # check dependencies
  library("purrr")
  library("data.table")
  
  # check if multithreading is required and make sure an integer of threads is used
  if(!is.null(p)) {
    if(!as.integer(p) == p) stop("Number of threads must be integer!")
    threads_old = setDTthreads(p)
  }
  
  # clean up the indir path and check whether it is valid
  indir = paste0(dirname(indir), "/", basename(indir), "/")
  if (!dir.exists(indir)) stop(paste("Directory not found:", indir))
  
  # obtain list of files in indir
  bed_files = list.files(path = indir, pattern = ".bed")
  
  if(! length(bed_files) > 0) stop("No files found.") # sanity check
  
  if(verbose == T) {
    cat("Reading in these files:\n")
    cat(bed_files, sep="\n")
    cat("\n")
  }
  
  # detect sample groups
  if (is.null(samples)) {
    samples = unique(sapply(strsplit(bed_files, sep), "[", 1))
  }
  
  if(verbose == T) {
    cat("Selected samples:\n")
    cat(samples, "\n")
    cat("\n")
  }
  
  samples = factor(samples) # convert to factor
  
  # detect replicates and convert to character vector
  replicates = lapply(as.character(samples), function(x) {
    bed_files %>%
      .[grepl(x, .)] %>% 
      sub(".bed", "", .) %>%
      sub(paste0(x,sep), "", .) %>% 
      paste0("rep",.)
  })
  
  names(replicates) = samples
  replicates = replicates[levels(samples)] # possibly reorder samples
  
  files  = sapply(as.character(samples), simplify = F, USE.NAMES = T,
                  function(x) list.files(path = indir, pattern = x))   # detect files
  names(files) = samples
  files = files[levels(samples)]                                       # make sure to adhere to right order
  files = lapply(files, function(x) x[grepl(".bed", x)])               # make sure to onlöy read bed files
  
  if(any(sapply(files, length)==0)) stop("At least one group is empty.") # sanity check
  
  # loop over files and read them in
  rawbeds = lapply(1:length(files), function(x) {
    
    lapply(1:length(files[[x]]), function(y) {
      
      if(verbose == T) cat("Reading in", files[[x]][[y]], "\n")
      
      
      # deal with different bed formats
      header_of_file = fread(paste0(indir, "/", files[[x]][[y]]), nrow = 1, showProgress = F)
      
      number_of_columns_in_file = ncol(header_of_file)
      
      # we are not interested in the names of reads and discard them
      # special case: only chromsome, start and count known (three columns).
      # in this case, we will use end = start+1 to simulate the read at 5'
      if(number_of_columns_in_file == 6) {
        columns_to_read = c(1:3,5)
      } else if (number_of_columns_in_file == 4) {
        columns_to_read = 1:4
      } else if (number_of_columns_in_file == 3) {
        columns_to_read = 1:3 
      }
      
      # actually read in file
      file_dt = fread(paste0(indir, "/", files[[x]][[y]]), select = columns_to_read, showProgress = verbose)
      
      if (number_of_columns_in_file == 3) { # we do not know the end position, so use 5'
        setnames(file_dt, c("chr", "start", "count"))
        file_dt[, end := start+1]
      } else {
        setnames(file_dt, c("chr", "start", "end", "count"))
      }
      
      
      setkey(file_dt, chr, start, end, count)
      
      # now convert our data.table to GRanges for clean working
      file_dt = with(file_dt, GRanges(seqnames = chr, IRanges(start = start, end = end), count = count))
      seqlevelsStyle(file_dt) = "UCSC" # deal with different styles being used
      file_dt = file_dt[order(seqnames(file_dt))]
      
      file_dt[]
    })
  })
  names(rawbeds) = names(files)
  
  # name replicates
  for (sample in names(rawbeds)) {
    names(rawbeds[[sample]]) = replicates[[sample]]
  }
  
  
  invisible(gc(verbose=F, full=F)) # clean up
  
  
  rawbeds
}

### Aim: helper function to detect whether we have multiple groups of samples present in
### a folder (think: multiple sequencing experiments)
### Example usage: detect_groups(indir)

detect_groups = function(indir, multiple_groups = F) {
  library("dplyr")
  
  ### sanity checks
  if(!dir.exists(indir)) stop("Invalid input directory.")
  if(!is.logical(multiple_groups)) stop("Invalid input for multiple_groups.")
  
  # are multiple groups expected? Then check which ones contain bed files
  
  # if multiple groups have not explicitly been set, do we have .bed files in indir?
  if( multiple_groups == F && (list.files(path = indir) %>% grepl(".bed", .) %>% any(.) )) {
    
    dirs = indir # if so, then simply return indir
    
  } else {       # otherwise check which directories contain bed files
    
    dirs = list.dirs(path = indir)
    
    contain_bed = sapply(1:length(dirs), function(x) { 
      list.files(path = dirs[x]) %>% grepl(".bed", .) %>% any(.) 
    })
    
    dirs = dirs[contain_bed] %>% paste0(., "/")
    
    if(length(dirs) == 0 | is.null(dirs)) stop("No valid directory detected.") # sanity check
    
  }
  
  dirs = sub("//", "/", dirs) # clean dir paths
  
  return(dirs)
}


### Aim: starting from a nested list created using load_bed_from_dir or _groups
### create annotated objects for each replicate, maintaining the same list structure
### Example usage: annotate_mapped_reads(indir)

annotate_mapped_reads = function(bed_list, annotation_file = NULL, annotation_object = NULL, add_missing = T,
                                 verbose = T, upstream = NULL, downstream = NULL,
                                 p = NULL, valid_chr = c(paste0("chr", 1:22), "chrX")) {
  # check dependencies
  library("dplyr")
  library("GenomicRanges")
  library("data.table")
  library("tools")
  library("rtracklayer")
  
  # check if multithreading is required and make sure an integer of threads is used
  if(!is.null(p)) {
    if(! as.integer(p) == p) stop("p must be integer")
  } else p = 1L
  
  # do we want to look at areas around the start/end of our reads?
  if(any(!is.null(upstream), !is.null(downstream)) ) {
    if(any(is.null(upstream), is.null(downstream))) stop("Please set both, upstream and downstream.")
    promoters = T 
  } else promoters = F
  
  # obtain annotation objection, either by reading it in or by taking it from arguments
  if(is.null(annotation_object)) {
    if(!is.null(annotation_file)) {
      if(verbose) cat("Reading in annotations ...")
      if (file_ext(annotation_file) == "bed") { all_regions = import.bed(annotation_file) 
      } else if (file_ext(annotation_file) == "rds") { all_regions = readRDS(annotation_file)
      } else stop("Unrecognized annotation file format. Accepted are bed and rds.")
    } else stop("Provide either annotation file or object.")
  } else {
    all_regions = annotation_object
  }
  
  # determine the column name of our count data
  count_col = unique(rapply(bed_list, function(x) names(x@elementMetadata)))
  if(! length(count_col) == 1) stop("Bed list is inconsistent with naming for count data!")
  
  # check whether we are looking a multiple groups, or just one
  l_depth = list_depth(bed_list)
  if(! ( l_depth == 2L | l_depth == 3L ) )  stop(paste("dt_list must either be depth 2 or 3, is", l_depth))
  
  # convert annotations to data.table format
  all_regions = as.data.table(all_regions)
  
  if(!all(c("seqnames", "start", "end") %in% colnames(all_regions))) { # sanity check
    stop("Annotations must contain seqnames, start and end columns!") }
  
  to_sort = colnames(all_regions)
  to_sort = to_sort[! to_sort %in% c("start", "end") ] # remove these two while collpasing
  
  # collapse annotation to remove duplicates and invalid chromosomes
  all_regions = all_regions[seqnames %in% valid_chr,
                            .(start = min(start), end = max(end)),
                            by = to_sort]
  
  to_sort = colnames(all_regions) # obtain column names for later use
  
  
  all_regions = makeGRangesFromDataFrame(all_regions, keep.extra.columns = T) # convert annotations back to GRanges
  
  
  if(promoters) { all_regions = promoters(all_regions, upstream = upstream, downstream = downstream) } # get promoters if wanted
  
  if(verbose) cat("Overlapping reads with annotation...")
  
  fo = rapply(bed_list, how = "list", function(x) findOverlaps(x, all_regions) ) # overlap reads with annotations
  
  if(l_depth == 2) { # groups do not matter, so do not introduce a separate column for them
    # bind together the sequence names and counts per gene
    dsb_in_regions = mclapply(mc.cores = p, names(fo),
                              function(x) {
                                lapply( names( fo[[x]] ), function(y) {
                                  
                                  # first, from the GRanges get the library size by the summation of all counts
                                  # this is done via accessing the metadata at the count column determined earlier
                                  lib_size = sum(bed_list[[x]][[y]]@elementMetadata[[count_col]])
                                  
                                  # second, create a data.table with the annotations(left) and our read info (right)
                                  tmp = cbind(
                                    as.data.table( all_regions[subjectHits(fo[[x]][[y]])] ) [,..to_sort],
                                    as.data.table( bed_list[[x]][[y]][queryHits(fo[[x]][[y]])] ) [,..count_col] 
                                  )
                                  
                                  to_sort_tmp = to_sort[! to_sort == "width" ] # we need width separate
                                  
                                  # third, calculate I) breaks per kb, II) counts per million (CPM) III) reads per kb and million reads (RPKM)
                                  tmp[, genadj := get(count_col)/width * 1000, by = to_sort_tmp]       # breaks per kb
                                  tmp[, width := max(width), by = to_sort_tmp] 
                                  tmp[, libadj := sum(get(count_col))/lib_size*1e06, by = to_sort_tmp] # CPM
                                  tmp[, libgenadj := sum(genadj)/lib_size*1e06, by = to_sort_tmp]      # RPKM
                                  tmp[, genadj := sum(genadj), by = to_sort_tmp]                       # collapse
                                  tmp[, eval(count_col) := sum(get(count_col)), by = to_sort_tmp]      # collapse
                                  
                                  # fourth, if wanted, add genes not covered in some samples
                                  if(add_missing) {
                                    to_add = cbind(
                                      as.data.table(all_regions)[! gene_id %in% unique(tmp$gene_id), ..to_sort],
                                      count = 0, genadj = 0, libadj = 0, libgenadj = 0
                                    )
                                    tmp = rbindlist(list(tmp, to_add), use.names = T)
                                  }
                                  
                                  setkeyv(tmp, to_sort_tmp) # sort the data.table
                                  
                                  tmp = unique(tmp)         # make sure to avoid duplicates
                                  setcolorder(tmp, c("seqnames", "start", "end", "width")) # priorise these columns
                                })
                              })
  } else { # groups do matter, so do introduce a separate column for them
    
    dsb_in_regions = mclapply(mc.cores = p, names(fo),
                              function(x) {
                                lapply( names( fo[[x]] ), function(y) {
                                  tmp_list = 
                                    lapply( names( fo[[x]][[y]] ), function(z) {
                                      
                                      # first, from the GRanges get the library size by the summation of all counts
                                      # this is done via accessing the metadata at the count column determined earlier
                                      lib_size = sum(bed_list[[x]][[y]][[z]]@elementMetadata[[count_col]])
                                      
                                      # second, create a data.table with the annotations(left) and our read info (right)
                                      tmp = cbind(
                                        as.data.table( all_regions[subjectHits(fo[[x]][[y]][[z]])] ) [,..to_sort],
                                        as.data.table( bed_list[[x]][[y]][[z]][queryHits(fo[[x]][[y]][[z]])] ) [,..count_col] 
                                      )
                                      
                                      to_sort_tmp = to_sort[! to_sort == "width" ] # we need width separate
                                      
                                      # third, calculate I) breaks per kb, II) counts per million (CPM) III) reads per kb and million reads (RPKM)
                                      tmp[, genadj := get(count_col)/width * 1000, by = to_sort_tmp]       # breaks per kb
                                      tmp[, width := max(width), by = to_sort_tmp] 
                                      tmp[, libadj := sum(get(count_col))/lib_size*1e06, by = to_sort_tmp] # CPM
                                      tmp[, libgenadj := sum(genadj)/lib_size*1e06, by = to_sort_tmp]      # RPKM
                                      tmp[, genadj := sum(genadj), by = to_sort_tmp]                       # collapse
                                      tmp[, eval(count_col) := sum(get(count_col)), by = to_sort_tmp]      # collapse
                                      
                                      # fourth, if wanted, add genes not covered in some samples
                                      if(add_missing) {
                                        to_add = cbind(
                                          as.data.table(all_regions)[! gene_id %in% unique(tmp$gene_id), ..to_sort],
                                          count = 0, genadj = 0, libadj = 0, libgenadj = 0
                                        )
                                        tmp = rbindlist(list(tmp, to_add), use.names = T)
                                      }
                                      
                                      setkeyv(tmp, to_sort_tmp) # sort the data.table
                                      
                                      tmp = unique(tmp) # make sure to avoid duplicates
                                      setcolorder(tmp, c("seqnames", "start", "end", "width")) # priorise these column
                                    })
                                  
                                  names(tmp_list) = names(bed_list[[x]][[y]]) # name innermost list
                                  tmp_list
                                  
                                })
                              })
  }
  
  # name list
  names(dsb_in_regions) = names(bed_list)
  
  for(sample in names(dsb_in_regions)) {
    names(dsb_in_regions[[sample]]) = names(bed_list[[sample]])
  }
  
  # cleanup RAM
  invisible(gc(full = F))
  
  dsb_in_regions
  
}


### Aim: wrapper for load_bed_from_dir for multiple groups. Additionally, offers possibility to annotate reads.
### Example usage: load_bed_from_groups(indir) [where indir contains files of pattern "NES_1.bed"]

load_bed_from_groups = function(indir, multiple_groups = F, annotation_file = NULL, annotation_object = NULL,
                                return_raw = F, upstream = NULL, downstream = NULL, sep = "_",
                                valid_chr = c(paste0("chr", 1:22), "chrX"), verbose = F, quiet = F, p = 1L) {
  # check dependencies
  library("data.table")
  library("parallel")
  
  
  if(quiet) verbose = F # quiet takes precedence over verbose
  
  
  ## sanity checks
  
  if(!dir.exists(indir)) stop("Invalid input directory.")
  
  if(!is.null(p)) {
    if(!as.integer(p) == p) stop("Number of threads must be integer!")
  } else p = 1
  
  # if annotation is wanted, determine if via object or file, and if needed provide a default
  if(!return_raw) {
    if(is.null(annotation_object)) {
      if(is.null(annotation_file) || !file.exists(annotation_file) ) {
        #annotation_file = "/media/jesko/HDD_Ubuntu/BLISS/annotations/hg19/general/hg19_Gencode19_annotations.all.genes.regions.rds"
        annotation_file = "/Users/gustaweriksson/Dokument/Bienko-Crosetto/BLISS_RNA_Pipeline/Data/hg19/annotation/hg19_Gencode19_annotations.all.genes.regions.rds"
        message(paste0("Defaulting annotation file to ", annotation_file,"."))
      }
    } else annotation_file = NULL # provided object takes precendence over file path
  }
  
  ## end sanity checks
  
  
  # determine which directories to look in for bed files
  groups = detect_groups(indir = indir, multiple_groups = multiple_groups)
  
  if(length(groups) == 1 && multiple_groups == T) warning("Multiple groups mode enabled, but finding only one group!")
  
  dsb_in_groups = lapply(groups, function(current_group) {
    if (!quiet) cat("Reading in", basename(current_group),"...\n")
    rawbeds = load_bed_from_dir(current_group, verbose = verbose, sep = sep, p = p)     # load in bed files
    samples = names(rawbeds)                                                            # detect samples
    rawbeds = rapply(rawbeds, how = "list", function(x) x[seqnames(x) %in% valid_chr])  # make sure to only have valid chromosomes
    
    if(return_raw) return(rawbeds)  # if only the nested list is required without annotation
    
    # otherwise overlap the bed files with gene annotations
    if(!quiet) cat("Annotating", basename(current_group),"...\n\n")
    annotated_sample = annotate_mapped_reads(bed_list = rawbeds,
                                             annotation_file = annotation_file,
                                             annotation_object = annotation_object,
                                             upstream = upstream, 
                                             downstream = downstream,
                                             p = p,
                                             verbose = verbose)
    
    return(annotated_sample)
  }) # lapply will merge all groups into a list
  
  names(dsb_in_groups) = basename(groups) # name list elements
  
  invisible(gc()) # clean up
  
  dsb_in_groups # return list
}



### Aim: After having loaded data from a directory using load_bed_from_dir/groups, and possibly annotating
### the reads using annotate_mapped_reads, make_DT_from_list converts the reads to an easy-to-handle data.table.
### This has the advantage of being much easier to work with for further processing.
### Example usage: make_DT_from_list(dt_list) [where dt_list contains a nested list of data.table or GRanges to be converted]

make_DT_from_list = function(dt_list) {
  
  # check dependencies
  library("data.table")
  library("GenomicRanges")
  library("purrr")
  
  ### sanity checks ###
  if(!is.list(dt_list)) stop("Provided object is not a list")
  
  
  l_depth = list_depth(dt_list)
  if(! ( l_depth == 2L | l_depth == 3L ) )  stop(paste("dt_list must either be depth 2 or 3, is", l_depth))
  
  if(all(rapply(dt_list, function(x) "GRanges" %in% class(x)) == T)) { 
    list_mode = "gr"
  } else {
    if(l_depth == 2L) {
      if("data.table" %in% class(dt_list[[1]][[1]])) {
        list_mode = "dt"
      }
    } else if(l_depth == 3L) {
      if("data.table" %in% class(dt_list[[1]][[1]][[1]])) { 
        list_mode = "dt"
      }
    }
  }
  
  if(!exists("list_mode")) stop("Could not recognize data type. Must be data.table or GRanges.")
  ### end sanity checks ###
  
  # go through the nested list, add meta-information about groups, samples, and replicate number,
  # then make one big data.table out of it and order that
  
  if(list_mode == "gr") {
    dt_list = rapply(dt_list, how = "list", f = as.data.table) # if using GRanges they will be converted
  }
  
  if(l_depth == 3L) {
    dt_list = 
      rbindlist(
        lapply(names(dt_list), function(group) {
          rbindlist(
            lapply(names(dt_list[[group]]), function(sample) {
              rbindlist(
                lapply(names(dt_list[[group]][[sample]]), function(replicate) {
                  dt_list[[group]][[sample]][[replicate]][, ':=' (group = group, sample = sample, replicate = replicate)]
                })
              )
            })
          )
        })
      )
  } else {
    dt_list = 
      rbindlist(
        lapply(names(dt_list), function(sample) {
          rbindlist(
            lapply(names(dt_list[[sample]]), function(replicate) {
              dt_list[[sample]][[replicate]][, ':=' (sample = sample, replicate = replicate)]
            })
          )
        })
      )
  }
  
  
  # specify order of columns to be returned
  
  if(l_depth == 3L) { 
    setcolorder(dt_list, c("group", "sample", "replicate"))
  } else {
    setcolorder(dt_list, c("sample", "replicate"))
  }
  
  
  # It can happen that in some samples, different "types" of the same gene id are found, 
  # resulting in a different width for the gene.
  # In these cases, one has to readjust the normalisation for gene length.
  # This happens in a bit less than 5% of all genes and only if the object you provided initially
  # also had the type column (which causes multiple entries for the same gene id).
  # Requiring further testing, this is currently not enabled by default.
  
  # if("width" %in% colnames(bed_list)) {
  #   tmp = bed_list[,.(n_width = length(unique(width))), by = gene_id][n_width > 1,gene_id]
  #   
  #   bed_list[gene_id %in% tmp, ':=' (width = max(width), libgenadj = libgenadj*width/max(width)), by = gene_id]
  # }
  
  invisible(gc())
  return(dt_list[]) # return
}


### Aim: Use genome binning to count total reads in a bin, or unique positions found in bins.
### This can later be used for e.g. correlation analysis.
### Example usage: bin_reads(bed_list, mode = "count") [where bed_list contains a nested list of GRanges to be used]

bin_reads = function(bed_list, binsize = NULL, mode = c("count", "unique"), samples = NULL, provided_bins = NULL, annotation_dir = NULL,
                     blacklist_file = NULL, adj_lib = T, collapse_rep = T, valid_chr = c(paste0("chr", 1:22), "chrX")) {
  
  # check dependencies
  library("data.table")
  library("GenomicRanges")
  library("rtracklayer")
  
  
  l_depth = list_depth(bed_list)
  if(! ( l_depth == 2L | l_depth == 3L ) )  stop(paste("bed_list must either be depth 2 or 3, is", l_depth))
  
  if(length(mode) == 0) {
    mode = "count"
    message("Defaulting mode to count.")
  } else if (length(mode) > 1) {
    mode = mode[1] # only the first choice is respected
    message(paste("Only first mode argument considered. mode =", mode))
  }
  
  if(!mode %in% c("count", "unique")) stop(paste("mode must be \"count\" or \"unique\", is"), mode)
  
  if(is.null(annotation_dir)) { annotation_dir  = "/media/jesko/HDD_Ubuntu/BLISS/annotations/hg19/general/"
  }
  if(!dir.exists(annotation_dir)) stop("Invalid annotation directory.")
  
  if(is.null(blacklist_file)) { blacklist_file = "/media/jesko/HDD_Ubuntu/BLISS/annotations/hg19/blacklist/broad_2k_merged.bed"
  }
  if(!file.exists(blacklist_file)) stop("Invalid blacklist file.")
  
  if(!is.null(samples)) samples = factor(samples, levels = samples) 
  
  if(!is.null(provided_bins)) {
    if("GRanges" %in% class(provided_bins)) {
      bins = provided_bins
    } else  { bins = makeGRangesFromDataFrame(provided_bins, keep.extra.columns = T) }
  } else if(binsize != as.integer(binsize)) stop("Binsize must be integer if no bins provided.") 
  
  
  # create bins from chromosme information, if no bins have been provided
  if(is.null(provided_bins)) {
    blacklist    = import.bed(blacklist_file) # use GenomicRanges to import for sanity
    chromosomes  = genomeStyles()[["Homo_sapiens"]]$UCSC # deal with various genome styles
    
    # get metainfo on chromosomes from files
    chrominfo    = sortSeqlevels(with(fread(paste0(annotation_dir, "hg19", ".chrom.sizes"))[V1%in%chromosomes,],
                                      Seqinfo(V1, seqlengths = V2, genome = "hg19")))
    chrsizes     = seqlengths(chrominfo)
    chrsizes     = chrsizes[names(chrsizes) %in% valid_chr]
    chrsizes     = GRanges(seqnames = names(chrsizes), ranges = IRanges(start=1, end=chrsizes))
    chrsizes_dt  = as.data.table(chrsizes)
    
    bins = lapply(1:nrow(chrsizes_dt), function(x) { # go along chromosomes and bin genome
      data.table(
        seqnames = chrsizes_dt[x,seqnames],
        start = seq(1, chrsizes_dt[x,end], by = binsize),
        end = c(seq(binsize, chrsizes_dt[x,end], by = binsize), chrsizes_dt[x,end])
      )
    }) %>% rbindlist %>% makeGRangesFromDataFrame() %>%
      subsetByOverlaps(., blacklist, invert = T) # remove any bin overlapping with blacklist
  }
  
  # go through list, assign their score (as read in) and find their overlap with the bins
  binned_reads = rapply(bed_list, how = "list", function(x) {
    if(mode == "count") {                                                                            # if we want all reads to be counted
      co    = countOverlaps(bins, x)                                                                 # count overlaps between reads and bins
      x     = data.table(chr=as.vector(seqnames(bins)), start=start(bins), end=end(bins), count=co)  # convert to data.table
    } else {                                                                                         # or if we only want unique reads
      fo    = findOverlaps(bins, x)                                                                  # then just determine overlaps
      x     = as.data.table(bins[queryHits(fo)])                                                     # convert to data.table
      x     = x[,.(count = .N), by=.(seqnames, start, end)]                                          # count unique positions
      if(nrow(x) != length(bins)) {                                                                  # add those bins where no read was found
        x   = rbind(x, as.data.table(bins[-queryHits(fo)])[,.(seqnames = seqnames, start = start, end = end, count = 0)])
      }
      setnames(x, "seqnames", "chr")
    }
  })
  
  # bind all the replicates together
  cts_binned = make_DT_from_list(binned_reads)
  cts_binned[, chr := factor(chr, levels = valid_chr)] # make sure to have proper chromosome factor levels
  setkey(cts_binned, sample, chr, start, end)
  
  if(adj_lib) {                                        # if wanted, the read count can be adjusted by library size
    if(!mode == "count") {
      warning("Adjusting by library size with mode other than count does not adjust for all reads but count of unique positions.")
    } 
    # adjust for library size (CPM)
    if(l_depth == 2L) { cts_binned[, libadj := count/(sum(count))*1e6, by = .(sample, replicate)]
    } else cts_binned[, libadj := count/(sum(count))*1e6, by = .(group, sample, replicate)]
  }
  
  # collapse replicates, if wanted, by taking the mean across replicates
  if(collapse_rep) {
    cts_binned[, replicate := NULL]
    cts_binned[, count := as.numeric(count)] # convert to numeric from integer to avoid warning
    if(l_depth == 2L) { cts_binned = cts_binned[, .(count = mean(count), libadj = mean(libadj)), by = key(cts_binned)]
    } else cts_binned = cts_binned[, .(count = mean(count), libadj = mean(libadj)), by = c("group", key(cts_binned))]
  }
  
  # add sample factor levels, if wanted
  if(!is.null(samples)) {
    cts_binned[, sample := factor(sample, levels = samples)]
  }
  
  setkeyv(cts_binned, colnames(cts_binned)) # sort from left to right
  
  return(cts_binned[])
}


### Aim: Caculate where breaks occur: in genic, or intergenic regions.
### This can later be used for plotting of pie charts with information on where breaks occur.
### Example usage: count_break_location(indir)

count_break_location = function(bed_list, annotation_file = NULL, groups = NULL, samples = NULL, 
                                valid_chr = c(paste0("chr", 1:22), "chrX"), quiet = F) {
  # check dependencies
  library("data.table")
  library("GenomicRanges")
  
  ### sanity checks and defaults
  if(is.null(annotation_file)) {
    annotation_file = "/media/jesko/HDD_Ubuntu/BLISS/annotations/hg19/general/hg19_Gencode19_annotations.all.genes.regions.rds"
    message(paste0("Defaulting annotation file to ", annotation_file,"."))
  }
  
  l_depth = list_depth(bed_list)
  if(! ( l_depth == 2L | l_depth == 3L ) )  stop(paste("bed_list must either be depth 2 or 3, is", l_depth))
  
  
  if(!file.exists(annotation_file)) stop("Annotation file not found")
  ### end sanity checks
  
  # read in annotations
  all_regions = readRDS(annotation_file)
  all_regions = as.data.table(all_regions)  # collapse all genes in annotation by gene_id
  all_regions[, ':=' (start = min(start), end = max(end)),
              by = .(seqnames, strand, gene_id, gene_name, annotation, type)]
  all_regions = all_regions[seqnames %in% valid_chr]
  all_regions = makeGRangesFromDataFrame(all_regions, keep.extra.columns = T) # make GRanges
  
  # clean chromosomes
  bed_list = rapply(bed_list, how = "list", function(x) x[seqnames(x) %in% valid_chr])
  
  # overlap
  ann_beds = rapply(bed_list, how = "list", function(x) subsetByOverlaps(x, all_regions))
  
  # count portions
  all_reads         = rapply(bed_list, how = "unlist", function(x) sum(x$count)) # determine raw counts
  reads_in_anno     = rapply(ann_beds, how = "unlist", function(x) sum(x$count))
  dif_count         = data.table(identifier = names(all_reads),
                                 all_reads = all_reads,
                                 reads_in_anno = reads_in_anno)
  
  # split the identifiers by the "." delimiter introduced by rapply
  # then calculate percentages
  if(l_depth == 2) {
    dif_count[, c("sample", "replicate") := tstrsplit(identifier, ".", fixed=T)]
    dif_count[,identifier := NULL] # remove unnecessary column
    dif_count[, ':=' (int_per = (all_reads-reads_in_anno)/(all_reads)*100, gen_per = reads_in_anno/(all_reads)*100)]
    setcolorder(dif_count, c("sample", "replicate", "all_reads", "reads_in_anno", "int_per", "gen_per"))
  } else {
    dif_count[, c("group", "sample", "replicate") := tstrsplit(identifier, ".", fixed=T)]
    dif_count[,identifier := NULL] # remove unnecessary column
    dif_count[, ':=' (int_per = (all_reads-reads_in_anno)/(all_reads)*100, gen_per = reads_in_anno/(all_reads)*100)]
    setcolorder(dif_count, c("group", "sample", "replicate", "all_reads", "reads_in_anno", "int_per", "gen_per"))
  }
  
  
  if(any(dif_count[,int_per+gen_per] != 100)) stop("Error in calculating percentages.") #sanity check 
  
  
  # set factors if wanted
  if(!is.null(groups)) {
    if(all(dif_count$sample %in% groups)) {
      dif_count[,group := factor(group, levels = groups)]
    } else {
      warning("Group detected that is not in the passed group vector.")
      print(dif_count[!group %in% groups, group])
      print(groups)
    }
  }
  
  
  if(!is.null(samples)) {
    if(all(dif_count$sample %in% samples)) {
      dif_count[,sample := factor(sample, levels = samples)]
    } else {
      warning("Group detected that is not in the passed sample vector.")
      print(dif_count[!sample %in% samples, sample])
      print(samples)
    }
  }
  
  if(l_depth == 2) {
    setkey(dif_count, sample, replicate)
  } else {
    setkey(dif_count, group, sample, replicate)
  }
  
  return(dif_count[])
}


### Aim: After obtaining a data.table holding information of bed files, one may wish to collapse
### reads from replicates.
### Examples usage: collapse_replicates_in_DT(dt, fun = sum())

collapse_replicates_in_DT = function(dt, fun = mean(), rep_colname = "replicate", exclude = c("start", "end", "width")) {
  # dependencies
  library("data.table")
  library("dplyr")
  
  # sanity checks
  if(!"data.frame" %in% class(dt)) stop("dt must be data.frame or data.table.")
  if(!rep_colname %in% colnames(dt)) stop(
    sprintf("Column %s not found within dt column names: %s", rep_colname, paste(colnames(dt), collapse = " "))
  )
  
  # check which columns to transform and which to group by
  numeric_columns = sapply(dt, is.numeric)
  count_cols = names(numeric_columns[numeric_columns == T]) %>% .[!. %in% exclude]
  group_cols = colnames(dt)[!colnames(dt) %in% count_cols] %>% .[!. == rep_colname]
  
  # sanity check
  if(is.null(count_cols) | is.null(group_cols)) {
    stop("No count or group column detected.")
  }
  
  # convert possible integer count columns to double to be able to take the mean
  integer_columns = sapply(dt, is.integer) 
  integer_columns = names(integer_columns)[integer_columns]
  integer_columns = integer_columns[!integer_columns %in% exclude]
  
  if(length(integer_columns) != 0) {
    dt[,eval(integer_columns) := as.double(get(integer_columns))]
  }
  
  return(dt[, lapply(.SD, mean), by=group_cols, .SDcols=count_cols])
}


### Aim: calculate the number of bases between two adjacent reads and return
### them as a data.table across samples and replicates
### Example usage: distances_between_reads(bed_list)

distances_between_reads = function(bed_list, valid_chr = c(paste0("chr", 1:22), "chrX")) {
  # dependencies
  library("GenomicRanges")
  library("data.table")
  
  # sanity check
  if(!is.list(bed_list)) stop("bed_list is not a list.")
  
  # recurse through list, select valid chromosomes, then get distances and make new GRanges
  to_return = rapply(bed_list, how = "list", function(x) {
    x = x[seqnames(x) %in% valid_chr]
    distances = distanceToNearest(x)
    GRanges(x[queryHits(distances)], distance = elementMetadata(distances, use.names = F))
  })
  
  to_return = make_DT_from_list(to_return) # convert to data.table for easier processing
  
  setnames(to_return, "distance.distance", "distance") # during creation of GRanges this name is created, might need fixing
  
  to_return[distance := distance+1] # ensure that adjacent reads have a distance of 1
  
  return(to_return[])
}


### function to obtain nice labels for ggplot2. Taken from: https://stackoverflow.com/a/56449202
### All credit to Luksurious.
number_format = function(x, digits = 0, space_delimited = T) {
  intl = c(1e3, 1e6, 1e9, 1e12)
  suffixes = c(' K', ' M', ' B', ' T')
  
  i = findInterval(x, intl)
  
  i_neg = findInterval(-x, intl)
  
  result = character(length(x))
  
  # Note: for ggplot2 the last label element of x is NA, so we need to handle it
  ind_format = !is.na(x) & i > 0
  neg_format = !is.na(x) & i_neg > 0
  
  # Format only the elements that need to be formatted
  # with suffixes and possible rounding
  result[ind_format] = paste0(
    formatC(x[ind_format] / intl[i[ind_format]], format = "f", digits = digits),
    suffixes[i[ind_format]]
  )
  # Format negative numbers
  result[neg_format] = paste0(
    formatC(x[neg_format] / intl[i_neg[neg_format]], format = "f", digits = digits),
    suffixes[i_neg[neg_format]]
  )
  
  # To the rest only apply rounding
  result[!ind_format & !neg_format] = as.character(
    formatC(x[!ind_format & !neg_format], format = "f", digits = digits)
  )
  
  if(!space_delimited) result = gsub(" ", "", result)
  return(result)
}

### Aim: Rank genes based on fragility, with options to use normalisation
### for reads/kb and filter to obtain only genes of desired size

rank_genes_by_fragility = function(dt, count_cols = "count", group_cols = "sample", 
                                   order = c(-1L, 1L), ties.method = "min",
                                   normalised = F, gene_length = NULL) {
  # dependencies
  library("data.table")
  
  # sanity check
  if(! (all(count_cols %in% colnames(dt)) & all(group_cols %in% colnames(dt))))
    
    if(!is.null(gene_length)) {
      if(! all(as.integer(gene_length) == gene_length))  {
        stop("gene_length parameters must be integer.")
      }
      if(! all(c("start", "end") %in% colnames(dt)) ) stop("Need start and end column for size selection.")
      # end sanity checks
      
      # select genes of desired length
      lower_length = min(gene_length)
      upper_length = max(gene_length)
      
      dt = dt[between(end-start, lower_length, upper_length)]
      
    }
  
  dt = copy(dt) # make sure not to change calling environment
  # note: if using huge data.tables, you can ommit this call 
  # to change the data.table in your environment to save memory
  
  # rank genes
  dt[, rank := frank(- get(count_cols), ties.method = ties.method), by = group_cols]
  
  return(dt[])
}
