#!/usr/bin/env Rscript
# env_bin_abundance_analysis.R
# Copyright Jackson M. Tsuji, 2019
# Auto-generates taxa plots of abundances of genome bins
# See required R packages below.

#####################################################
## Load required packages: ##########################
library(getopt)
library(futile.logger)
library(plyr)
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(glue))
library(scales)
#####################################################

# #####################################################
# ## User variables - only used if not running from the console
# # See proper descriptions of variables under parse_command_line_input
# # To run, uncomment and then run main(opt)
# opt <- list()
# 
# # Required inputs
# # setwd("/home/jmtsuji/Research_General/PhD/04b_Metagenome_resequencing_F2015/10_ATLAS_re_analysis/09_env_abundance_of_bins/vs2_dRep_scripted")
# opt$input_bin_mapping_stats_filename <- "genome_bin_mapping_stats.tsv"
# opt$input_gtdb_classification_bac120_filename <- "ELA111314_dRep_gtdbtk.bac120.summary.tsv"
# opt$input_gtdb_classification_ar122_filename <- "ELA111314_dRep_gtdbtk.ar122.summary.tsv"
# opt$output_prefix <- "test1"
# 
# # Optional inputs (set to 'NA' to ignore)
# opt$taxonomic_rank <- "family"
# opt$abundance_threshold <- 10 # If >= 1, shows the top x most abundant groups. 
#                                 # If < 1, shows all group over x proportional abundance
# opt$plot_legend_taxonomy <- "all"
# 
# # Run the script AFTER sourcing
# # main(opt)
# #####################################################

# Function: loads and joins the stats and classification tables
load_and_join_tables <- function(input_bin_mapping_stats_filename, input_gtdb_classification_bac120_filename,
                                 input_gtdb_classification_ar122_filename) {
  
  bin_mapping_stats <- read.table(input_bin_mapping_stats_filename, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  gtdb_classifications_bac120 <- read.table(input_gtdb_classification_bac120_filename, sep = "\t", header = TRUE,
                                            stringsAsFactors = FALSE)
  
  if (is.na(input_gtdb_classification_ar122_filename) == FALSE) {
    gtdb_classifications_ar122 <- read.table(input_gtdb_classification_ar122_filename, sep = "\t", header = TRUE,
                                             stringsAsFactors = FALSE)
    
    ### Combine the bac120 and ar122, first checking there are no duplicates
    # Check if any of the entries are shared (each shared entry will show as a vector of length nrow(gtdb_classifications_bac120))
    shared_entry_check <- gtdb_classifications_bac120$user_genome %in% gtdb_classifications_ar122$user_genome
    if ((TRUE %in% shared_entry_check) == TRUE) {
      
      # How many are TRUE? What are some examples?
      shared_entries <- gtdb_classifications_bac120$user_genome[shared_entry_check %in% TRUE]
      num_shared_entries <- length(shared_entries)
      
      # Error message
      stop("The following ", num_shared_entries, " entries are shared between the two GTDB files: '", 
           glue::glue_collapse(shared_entries, sep = ", "), "'. Cannot continue.")
    }
    
    # Combine
    gtdb_classifications <- dplyr::bind_rows(gtdb_classifications_bac120, gtdb_classifications_ar122)
    
  } else {
    gtdb_classifications <- gtdb_classifications_bac120
  }
  
  # Change name of genome column to prepare for left_join
  gtdb_classifications <- dplyr::rename(gtdb_classifications, "genome" = user_genome)
  
  # Simplify GTDB table (reduce total columns)
  gtdb_classifications <- dplyr::select(gtdb_classifications, genome, classification)
  ## TODO - consider checking out the other entries in more detail
  
  ### Join to the bin mapping stats
  final_stats_table <- dplyr::left_join(bin_mapping_stats, gtdb_classifications, by = "genome")
  
  return(final_stats_table)
}


# Function: calculates some more meaningful bin abundance stats from the raw data
# Input: 'final_stats_table' - tibble generated from 'read_and_join_tables()'
# Return: 'final_stats_table' with new columns for the final stats, and raw stats removed. Same number of rows.
calculate_mapping_stats <- function(final_stats_table) {
  
  # Reads per kilobase genome per billion total metagenomic reads (RPKB; modified form of RPKM used for genomics)
  final_stats_table$rpkb <- final_stats_table$mapped_reads / final_stats_table$genome_length_bp * 1000 / 
    final_stats_table$metagenome_total_reads * 10^9
  
  final_stats_table$percent_read_recruitment <- final_stats_table$mapped_reads / 
    final_stats_table$metagenome_total_reads * 100
  
  # Select down to the most meaningful stats
  final_stats_table <- dplyr::select(final_stats_table, genome, metagenome, coverage_mean, coverage_sd, 
                                     rpkb, mapped_reads, metagenome_total_reads, percent_read_recruitment, 
                                     classification)
  
  return(final_stats_table)
  
}

# Function: Recursively identify ambiguity hits (i.e., occurances of TRUE) starting from the back of the taxonomy_ranks
# Inputs: 'ambiguous_ranks' - logical (length 7) of whether or not each position in a taxonomy_ranks vector is considered ambiguous
# 'rank_number' - numeric (length 1). Must start at the end of the ambiguous_ranks, and then the function will propogate back to the start
# Return: 'ambiguous_rank_position' - numeric (length 1) of the 'last' ambiguous position in the ranks when moving back from the end of the ranks
check_ambiguity_pattern <- function(ambiguous_ranks, rank_number) {
  if ( ambiguous_ranks[rank_number] == TRUE ) {
    if ( rank_number == 1 ) {
      # End point 3: this means that you've reached all the way to the front of the taxonomy ranks and all have been ambiguous.
      # Return "1" to the user
      ambiguous_rank_position <- rank_number
      return(ambiguous_rank_position)
      
    } else if ( rank_number > 1 && rank_number <= 7 ) {
      # Recursion: If you're still at the far right or the middle of the taxonomic ranks, then all is well. Keep moving back until you hit the end.
      # message("Found ambiguity at position ", rank_number, ". Moving back one.")
      check_ambiguity_pattern(ambiguous_ranks, (rank_number - 1))
      
    } else {
      stop("rank_number must be a number <= 7. You provided '", rank_number, "'. Exiting...")
    }
    
  } else if ( ambiguous_ranks[rank_number] == FALSE ) {
    
    if ( rank_number == length(ambiguous_ranks) ) {
      # End point 1: this means you've just started running from the end of the ranks and found no ambiguities on the far right side
      # Thus, no ambiguous ranks of the proper pattern are here.
      ambiguous_rank_position <- NA
      return(ambiguous_rank_position)
      
    } else if ( rank_number < length(ambiguous_ranks) ) {
      # End point 2: you've propagated up the taxonomy ranks and eventually hit a rank that is non-ambiguous. Flag this.
      ambiguous_rank_position <- rank_number + 1
      return(ambiguous_rank_position)
      
    } else {
      stop("rank_number must be a number <= 7. You provided '", rank_number, "'. Exiting...")
    }
    
  }
}


# Function: parses GTDB taxonomy
# Input: a length-1 character vector of a GTDB taxonomy entry
# Return: a length-7 character vector with each GTDB taxonomic rank as an entry
parse_GTDB_taxonomy <- function(taxonomy_entry, resolve = FALSE) {
  # EXAMPLE (note - no "kindgom"):
  # d__Bacteria;p__Bacteroidota;c__Chlorobia;o__Chlorobiales;f__Chlorobiaceae;g__Chlorobium;s__
  
  # Check input is appropriate
  if (length(taxonomy_entry) != 1 | is.vector(taxonomy_entry) == FALSE) {
    flog.error("The provided taxonomy entry is inappropriate -- should be a vector of length 1. Exiting...")
    exit(1)
  }
  
  # Separate by semicolon
  parsed_taxonomy <- strsplit(taxonomy_entry, ";")[[1]]
  
  # Remove info header (e.g., "d__")
  parsed_taxonomy <- unlist(lapply(parsed_taxonomy, 
                                   function(x) { substr(x, start = 4, stop = nchar(x)) }))
  
  if (resolve == FALSE) {
    parsed_taxonomy <- gsub(pattern = "^$", x = parsed_taxonomy, replacement = "Unresolved")
  } else if (resolve == TRUE) {
    
    # Find all occurances of the ambiguous_terms among the taxonomy_ranks (output is TRUE/FALSE for each taxonomy rank)
    ambiguous_ranks <- grepl(pattern = "^$", x = parsed_taxonomy)
    
    # Get the position in the taxonomic ranks where the chain of ambiguous names ends (starting from the back)
    ambiguous_rank_position <- check_ambiguity_pattern(ambiguous_ranks, rank_number = 7)
    
    # Check for orphaned ambiguous ranks (i.e., ones that don't propagate in a chain from the back)
    # Don't modify orphaned entries, but throw a warning
    first_ambiguous_position <- match(TRUE, ambiguous_ranks)
    if ( is.na(first_ambiguous_position) == TRUE ) {
      # Don't do anything -- there are no ambiguous ranks for this entry.
    } else if ( first_ambiguous_position < ambiguous_rank_position ) {
      
      warning("Ambiguous taxonomy terms don't occur in a natural chain moving up from 'species' for this taxonomy entry: '",
              glue::glue_collapse(parsed_taxonomy, "; "), "'. Will only correct for the entries in a chain at the end.")
      
    }
    
    # Append taxonomic information to ambiguous entries
    if ( is.na(ambiguous_rank_position) == TRUE ) {
      # All done -- no ambiguous entries
    } else if ( ambiguous_rank_position > 1 && ambiguous_rank_position <=7 ) {
      
      # Get the taxon rank immediately above the ambiguous_rank_position
      last_informative_taxon <- parsed_taxonomy[ambiguous_rank_position - 1]
      
      # Append to end of ambiguous ranks
      for ( i in ambiguous_rank_position:7 ) {
        parsed_taxonomy[i] <- paste("Unresolved_", last_informative_taxon, sep = "")
      }
    }
  }
  
  # Check length looks okay
  if (length(parsed_taxonomy) != 7) {
    flog.error("The provided taxonomy entry does not contain seven ranks: '", glue_collapse(parsed_taxonomy, sep = "; "),
               "'. Exiting...")
    exit(1)
  }
  
  return(parsed_taxonomy)
  
}


# Function: Helper function to quickly get the number of unique entries in a column of a tibble
count_unique_ranks <- function(stats_table, grouping_var) {
  number_of_ranks <- length(unique(dplyr::pull(ungroup(stats_table), grouping_var)))
  
  return(number_of_ranks)
}


# Function: collapses the stats table by taxonomic rank
# Inputs: 'stats_table_with_ranks' - tibble of 'final_stats_table' with parsed GTDB ranks; 
#         'rank_to_summarize' - one of c("domain", "phylum", "class", "order", "family", "genus", "species"), 
#                               OR NA to skip this step and summarize at bin level (will just output the input)
# Return: collapsed input table
summarize_to_taxonomy <- function(stats_table_with_ranks, rank_to_summarize = NA) {
  # HARD-CODED
  tax_ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  
  if( is.na(rank_to_summarize) == TRUE) {
    
    flog.info("Not summarizing by taxonomic rank. Returning results at bin-level.")
    return(stats_table_with_ranks)
    
  } else {
    
    # Check input is okay
    if ((rank_to_summarize %in% tax_ranks) == FALSE) {
      stop("You provided the taxonomic rank of '", rank_to_summarize, "' to summarize to, but this rank does not exist. ", 
           "Please provide a canonical taxonomic rank (i.e., one of '", glue::glue_collapse(tax_ranks, sep = ", "), "') ",
           "for summarizing, or 'NA' to skip summarizing. Exiting...")
    }
    
    flog.info(paste("Summarizing table to rank '", rank_to_summarize, "'", sep = ""))
    
    # Group the data by the desired taxonomic rank
    group_by_vars <- c("metagenome", tax_ranks[1:match(x = rank_to_summarize, table = tax_ranks)])
    stats_table_grp <- dplyr::group_by_at(stats_table_with_ranks, group_by_vars)
    
    # Sum the stats -- unfortunately, standard deviation is lost here
    stats_table_summ <- dplyr::summarise(stats_table_grp, coverage_mean = sum(coverage_mean, na.rm = TRUE),
                                         rpkb = sum(rpkb, na.rm = TRUE), percent_read_recruitment = 
                                           sum(percent_read_recruitment, na.rm = TRUE))
    
    # Calculate stats
    pre_collapse <- count_unique_ranks(stats_table_with_ranks, "genome")
    post_collapse <- count_unique_ranks(stats_table_summ, rank_to_summarize)
    
    flog.info(paste("Collapsed from", pre_collapse, "genome bins to", post_collapse, "taxonomic groups"))
    
    return(stats_table_summ)
    
  }
  
}


# Function: drops groups in the data below a given abundance threshold
# Input 'stats_table_summarized' - tibble of 'final_stats_table' collapsed to desired taxonomic rank, with parsed GTDB ranks
#       'abundance_threshold' - numeric value. If >= 1, shows the top x most abundant groups. If < 1, shows all group over x proportional abundance.
# Return: input tibble with some dropped rows
# TODO - let user choose the stat
drop_low_abundance_groups <- function(stats_table_summarized, abundance_threshold, rank_to_summarize) {
  
  # If user has not collapsed to a taxonomic rank, then summarize by genome bin name
  if (is.na(rank_to_summarize)) {
    rank_to_summarize <- "genome"
  }
  
  # Calculate the starting number of groups for comparison later
  pre_drop <- count_unique_ranks(stats_table_summarized, rank_to_summarize)
  
  if (abundance_threshold >= 1) {
    
    flog.info(paste("Keeping the top", abundance_threshold, "most abundant groups in the data"))
    stats_table_summarized <- dplyr::group_by(ungroup(stats_table_summarized), metagenome)
    stats_table_summarized <- dplyr::top_n(stats_table_summarized, n = abundance_threshold, wt = rpkb)
    
    post_drop <- count_unique_ranks(stats_table_summarized, rank_to_summarize)
    flog.info(paste("Reduced total number of groups from", pre_drop, "to", post_drop))
    
  } else if (abundance_threshold < 1) {
    
    flog.info(paste("Keeping all groups in the data above", abundance_threshold, "proportional abundance"))
    
    # TODO - this is a bit janky.
    # Rough percent relative abundance for rpkb by summing rpkb (does NOT consider unassembled reads!!!)
    rpkb_totals <- dplyr::group_by(stats_table_summarized, metagenome)
    rpkb_totals <- dplyr::summarise(rpkb_totals, rpkb_sum = sum(rpkb, na.rm = TRUE))
    stats_table_summarized <- dplyr::left_join(stats_table_summarized, rpkb_totals, by = "metagenome")
    stats_table_summarized$rpkb_perc <- stats_table_summarized$rpkb / stats_table_summarized$rpkb_sum
    stats_table_summarized$rpkb_sum <- NULL
    
    stats_table_summarized <- dplyr::filter(stats_table_summarized, rpkb_perc >= abundance_threshold)
    stats_table_summarized$rpkb_perc <- NULL
    
    post_drop <- count_unique_ranks(stats_table_summarized, rank_to_summarize)
    flog.info(paste("Reduced total number of groups from", pre_drop, "to", post_drop))
    
  } else if (is.na(abundance_threshold) == TRUE) {
    
    flog.info(paste("Not dropping any low abundance groups. Total groups remains at", pre_drop))
    
  }
  
  return(stats_table_summarized)
  
}


# Function: chooses which discrete colour scale would work best for a plot based on the number of entries to be plotted
choose_discrete_colour_scale <- function(number_of_entries) {
  
  # Choose the best colour scale based on the number of entries to be plotted
  if ( number_of_entries == 2 ) {
    colour_palette <- RColorBrewer::brewer.pal(n = 3, name = "Dark2")[c(1,2)]
  } else if ( number_of_entries <= 8 ) {
    colour_palette <- RColorBrewer::brewer.pal(n = number_of_entries, name = "Dark2")
  } else if ( number_of_entries <= 12 ) {
    colour_palette <- RColorBrewer::brewer.pal(n = number_of_entries, name = "Set3")
  } else if ( number_of_entries > 12 ) {
    colour_palette <- hue_pal(h = c(20,290))(number_of_entries)
  } else {
    flog.error(paste("Something is wrong with the number_of_entries ('", number_of_entries, "'). Is it non-numeric? Exiting...", sep = ""))
  }

  return(colour_palette)
  
}


#
# Function: generates a basic taxaplot
# Inputs: 'stats_table' - tibble in the general form of 'final_stats_table' or derivative.
# Return: basic ggplot (can overlay custom colours later)
make_base_barplot <- function(stats_table, taxonomic_rank, add_taxonomy, plot_y_axis = "rpkb") {
  # HARD-CODED
  tax_ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  
  # Check add_taxonomy option
  if (add_taxonomy %in% c("all", "genus", NA) == FALSE) {
    flog.error(paste("add_taxonomy should be either 'all', 'genus', or 'NA', but you provided '", add_taxonomy,
                     '. Exiting...'))
    exit(1)
  }
  
  # If user has not collapsed to a taxonomic rank, then plot by genome bin name
  if (is.na(taxonomic_rank) | taxonomic_rank == "genome") {
    
    # Make taxonomic rank name match the column name for genomes
    taxonomic_rank = "genome"
    
    # Map taxon info onto the genome name
    if (is.na(add_taxonomy) == TRUE) {
    } else if(add_taxonomy == "all") {
      stats_table$genome <- paste(stats_table$genome, " (d__", stats_table$domain, ";p__", stats_table$phylum, 
                                  ";c__", stats_table$class, ";o__", stats_table$order, ";f__", stats_table$family, 
                                  ";g__", stats_table$genus, ";s__", stats_table$species, ")", sep = "")
    } else if (add_taxonomy == "genus") {
      stats_table$genome <- paste(stats_table$genome, " (g__", stats_table$genus, ")", sep = "")
    }
    
    # Arrange groups alphabetically by hierachical taxonomic ordering down to species level
    ranks_used <- tax_ranks
    stats_table <- dplyr::group_by_at(ungroup(stats_table), ranks_used)
    stats_table <- dplyr::arrange(stats_table, .by_group = TRUE)
    factor_ordering <- unique(dplyr::pull(stats_table, taxonomic_rank))
    stats_table[taxonomic_rank] <- factor(dplyr::pull(stats_table, taxonomic_rank), 
                                          levels = factor_ordering, ordered = TRUE)
    
  } else if (taxonomic_rank %in% tax_ranks) {
    
    # Map taxon info onto the rank name
    if (is.na(add_taxonomy) == TRUE) {
    } else if(add_taxonomy == "all") {
      # TODO - make this more comprehensive
      stats_table[taxonomic_rank] <- paste(stats_table[taxonomic_rank], 
                                           glue_collapse(dplyr::select(stats_table, 
                                                      tax_ranks[1:match(taxonomic_rank, table = tax_ranks)], sep = ";")))
    } else if (add_taxonomy == "genus") {
      flog.error("Cannot use option 'genus' when plotting to a taxon rank. Exiting...")
      exit(1)
    }
    
    # Arrange groups alphabetically by hierachical taxonomic ordering
    ranks_used <- tax_ranks[1:match(taxonomic_rank, table = tax_ranks)]
    stats_table <- dplyr::group_by_at(ungroup(stats_table), ranks_used)
    stats_table <- dplyr::arrange(stats_table, .by_group = TRUE)
    factor_ordering <- unique(dplyr::pull(stats_table, taxonomic_rank))
    stats_table[taxonomic_rank] <- factor(dplyr::pull(stats_table, taxonomic_rank), 
                                          levels = factor_ordering, ordered = TRUE)
    
  } else {
    flog.error("Provided taxonomic rank is not one of the possible options. Exiting...")
    exit(1)
  }
  
  # Determine the best colour scale
  number_of_ranks <- nrow(unique(dplyr::select(ungroup(stats_table), taxonomic_rank)))
  flog.info(paste("Plotting", number_of_ranks, "unqiue groups"))
  colour_scale <- choose_discrete_colour_scale(number_of_ranks)
  
  # Make the plot
  taxaplot <- ggplot(stats_table, aes(x = metagenome)) +
    geom_bar(aes_string(weight = plot_y_axis, fill = taxonomic_rank)) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_fill_manual(values = colour_scale) +
    xlab("Sample") +
    ylab(paste(plot_y_axis))
  
  return(taxaplot)
  
}


# Summary function for loading the tables
read_GTDB_abundance_data <- function(input_bin_mapping_stats_filename, input_gtdb_classification_bac120_filename,
                                     input_gtdb_classification_ar122_filename = NA,
                                     clarify_unresolved_taxa = FALSE) {
  
  # Load the data
  flog.info("Loading data frames")
  stats_table <- load_and_join_tables(input_bin_mapping_stats_filename, 
                                            input_gtdb_classification_bac120_filename,
                                            input_gtdb_classification_ar122_filename)
  
  flog.info("Calculating statistics")
  stats_table <- calculate_mapping_stats(stats_table)
  
  # Run the taxonomy parser on all rows and then convert the list output to a data frame (a bit messy)
  flog.info("Parsing GTDB taxonomy")
  taxonomy <- as.data.frame(matrix(unlist(
    lapply(stats_table$classification, function(x) { parse_GTDB_taxonomy(x, resolve = clarify_unresolved_taxa) })), 
    nrow = nrow(stats_table), byrow = TRUE), stringsAsFactors = FALSE)
  colnames(taxonomy) <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  
  # Bind onto table
  stats_table <- dplyr::bind_cols(stats_table, taxonomy)
  stats_table$classification <- NULL
  
  return(stats_table)
  
}


# Summary function for summarizing stats
summarize_GTDB_stats <- function(stats_table, rank_to_summarize, abundance_threshold, plot_legend_taxonomy,
                                 plot_y_axis = "rpkb") {
  
  # Summarize to desired taxonomic rank
  stats_table_summarized <- summarize_to_taxonomy(stats_table, rank_to_summarize)
  
  # Omit low-abundance groups if desired
  stats_table_summarized <- drop_low_abundance_groups(stats_table_summarized, abundance_threshold, rank_to_summarize)
  
  # Make the plot
  generic_taxaplot <- make_base_barplot(stats_table_summarized, rank_to_summarize, plot_legend_taxonomy,
                                        plot_y_axis = plot_y_axis)
  print(generic_taxaplot)

  # Return to user
  return_list <- list(stats_table_summarized, generic_taxaplot)
  names(return_list) <- c("stats_table_summarized", "generic_taxaplot")
  
  return(return_list)
  
}


main <- function(opt) {
  
  # Startup messages
  flog.info("Running env_bin_abundance_analysis.R")
  flog.info(paste("Input bin mapping stats filepath:", opt$input_bin_mapping_stats_filename))
  flog.info(paste("Input GTDB classification filename (bac120):", opt$input_gtdb_classification_bac120_filename))
  flog.info(paste("Input GTDB classification filename (ar122):", opt$input_gtdb_classification_ar122_filename))
  flog.info(paste("Output prefix:", opt$output_prefix))
  flog.info(paste("Taxonomic tank to summarize to for plotting (ignored if 'NA'):", opt$taxonomic_rank))
  flog.info(paste("Abundance threshold for plotting (ignored if 'NA'):", opt$abundance_threshold))
  flog.info(paste("Whether to add group taxonomy to each legend entry ('all', 'genus', or 'NA'):", opt$plot_legend_taxonomy))
  
  GTDB_stats_data <- read_GTDB_abundance_data(opt$input_bin_mapping_stats_filename, 
                                              opt$input_gtdb_classification_bac120_filename,
                                              opt$input_gtdb_classification_ar122_filename)
  
  
  GTDB_stats_summarized <- summarize_GTDB_stats(GTDB_stats_data, opt$taxonomic_rank, opt$abundance_threshold,
                                                opt$plot_legend_taxonomy)
  
}


# Assign command line input to variables, or throw help message, if R is being run from the console
if (interactive() == FALSE) {
  
  # Define flags
  params <- matrix(c('tree_filepath', 'i', 1, "character",
                     'blast_table_filepath', 'j', 1, "character",
                     'output_filepath', 'o', 1, "character",
                     'tree_metadata_filename', 'm', 2, "character",
                     'tree_decorator_colname', 'd', 2, "character",
                     'genome_plotting_names', 'n', 2, "logical",
                     'gene_naming_table_filename', 'g', 2, "character",
                     'bootstrap_cutoff', 'b', 2, "character",
                     'root_name', 'r', 2, "character",
                     'help', 'h', 2, "character"), byrow=TRUE, ncol=4)
  
  opt <- getopt(params)
  
  # If help was called, print help message and exit
  if ( !is.null(opt$help) ) {
    
    cat("generate_BackBLAST_heatmap.R: Binds a phylogenetic tree to a BLAST table heatmap.\n")
    cat("Copyright Lee Bergstrand and Jackson M. Tsuji, 2018\n\n")
    
    cat(getopt(params, usage = TRUE))
    
    cat("\n")
    
    message(glue::glue("
                       Required inputs:
                       --tree_filepath               Filepath for newick-format phylogenetic tree of the BLAST subject organisms
                       --blast_table_filepath        Filepath for CSV-format BLAST hit table from CombineBlastTables.R
                       --output_filepath             Output filepath for the PDF
                       
                       Optional inputs:
                       --tree_metadata_filename      Filepath for TSV-format metadata file for the phylogenetic tree. Details below.
                       --tree_decorator_colname      Column name from the metadata to map onto the tree as fill colours.
                       Requires that --tree_metadata_filename is set. Details below.
                       --genome_plotting_names       Set this flag to include a vector of genome names to plot in the tree_metadata_filename
                       called 'plotting_name'. Details below.
                       --gene_naming_table_filename  Filepath for TSV-format table linking query gene IDs in the BLAST table to proper
                       gene names. Details below.
                       --bootstrap_cutoff            A percentage at or above which to display the bootstrap values on the tree (e.g., 80)
                       --root_name                   Exact name of the tip you wish to use as the root of the tree, if your tree is not 
                       already rooted
                       
                       Tree vs. the BLAST table
                       Note that the subject organism names MUST be EXACTLY the same between the tree tips and the BLAST table 
                       'subject_name' column. Otherwise, the script will fail.
                       
                       Tree metadata:
                       You can optionally provide a tree_metadata table to overlay additional information onto the phylogenetic tree.
                       
                       The first column of this TSV (tab-separated) table MUST be called 'subject_name' and include the EXACT names
                       of all genomes in the phylogenetic tree. You can then optionally include the following:
                       
                       1. genome_plotting_names: add a column called 'plotting_name' to include the names of the organisms that you 
                       want to appear on the final plot. You can use spaces, most special characters, and so on.
                       
                       2. tree_decorator_colname: you can overlay characteristics of the organisms/genomes as fill colours on the tips
                       of the tree. Add a column to the metadata table with any name you'd like, e.g., 'GC_content', 
                       'predicted_metabolism', 'pathogenicity', and so on. Fill with meaningful data to you. Then, specify the
                       'tree_decorator_colname' to EXACTLY match the name of ONE of the additional columns. That info will then
                       be plotted on the tree. Enjoy! (We might expand this in the future to allow for font colours and so on to 
                       be varied.)
                       
                       Gene naming for BLAST table:
                       You can provide a gene_naming_table to provide custom names and ordering for query genes in your BLAST search, 
                       in place of the qseqid from NCBI, which may not be very human-readable.
                       
                       The gene naming table must meet the following criteria:
                       - First column: 'qseqid' - the EXACT ID of ALL of the unique query proteins in the BLAST table must be included.
                       You can include additional qseqid's here if you'd like, they just won't be used.
                       - Second column: 'gene_name' - a corresponding name of your choice (e.g., rpoB, dsrA, and so on)
                       
                       The order of the rows in this table will dictate the order of the genes in the heatmap.
                       
                       "))
    
    quit(status = 1)
  }
  
  # Exit if required inputs are not provided
  if ( is.null(opt$tree_filepath) ) {
    stop("Filepath to phylogenetic tree required (-i). Try -h for help message.")
  } else if ( is.null(opt$blast_table_filepath) ) {
    stop("Filepath to BLAST hit table required (-j). Try -h for help message.")
  } else if ( is.null(opt$output_filepath) ) {
    stop("Output filepath required (-o). Try -h for help message.")
  }
  
  # Check on the tree_metadata_filename and tree_decorator_colname were provided, as needed
  if ( is.null(opt$tree_metadata_filename) ) {
    # Annul everything
    # TODO - throw a warning if some of the other flags were set
    opt$tree_metadata_filename <- NA
    opt$tree_decorator_colname <- NA
    opt$genome_plotting_names <- NA
  } else if ( is.null(opt$tree_metadata_filename) == FALSE && is.null(opt$genome_plotting_names) == FALSE
              && is.null(opt$tree_decorator_colname) == FALSE ) {
    # All are present
    opt$genome_plotting_names <- TRUE
  } else if ( is.null(opt$tree_metadata_filename) == FALSE && is.null(opt$genome_plotting_names) == FALSE ) {
    # tree_decorator_colname must not be set
    opt$tree_decorator_colname <- NA
    opt$genome_plotting_names <- TRUE
  } else if ( is.null(opt$tree_metadata_filename) == FALSE && is.null(opt$tree_decorator_colname) == FALSE ) {
    # genome_plotting_names must not be set
    opt$genome_plotting_names <- NA
  }
  
  # Set defaults
  if ( is.null(opt$gene_naming_table_filename) == TRUE ) {
    opt$gene_naming_table_filename <- NA
  }
  if ( is.null(opt$bootstrap_cutoff) == TRUE ) {
    opt$bootstrap_cutoff <- NA
  }
  if ( is.null(opt$root_name) == TRUE ) {
    opt$root_name <- NA
  }
  
  # Run the script
  main(opt)
  
}
