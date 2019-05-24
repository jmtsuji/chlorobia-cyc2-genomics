# Figure_03_plotter.R
# Generates Figure 3 for the Chlorobi cyc2 paper -- environmental abundance of MAGs
# Copyright Jackson M. Tsuji, 2018 (Neufeld Research Group)
# See required packages in the library header here as well as in the headers of the two dependency scripts below

##### Load libraries and dependency scripts
library(here)
library(futile.logger)
library(glue, warn.conflicts = FALSE)
library(tibble, warn.conflicts = FALSE)
library(plyr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(reshape2, warn.conflicts = FALSE)
library(ggplot2, quietly = TRUE, warn.conflicts = FALSE)
library(gridExtra, warn.conflicts = FALSE)
library(egg, warn.conflicts = FALSE)
#####

##########################
### User-defined variables
##########################
params <- list()

# GTDB and read mapping
params$gtdb_filepath_ar122 <- here::here("input_data", "ELA111314_dRep_gtdbtk.ar122.summary.tsv")
params$gtdb_filepath_bac120 <- here::here("input_data", "ELA111314_dRep_gtdbtk.bac120.summary.tsv") 
params$mapping_stats_filename <- here::here("input_data", "genome_bin_mapping_stats.tsv")
params$metagenome_naming_info_filename <- here::here("input_data", "metagenome_naming_info.tsv")
params$assembled_read_stats_filename <- here::here("input_data", "assembled_read_stats.tsv")
params$checkM_stats_filepath <- here::here("input_data", "checkm_stats_reduced.tsv")
params$clarify_unresolved_taxa <- FALSE

# MetAnnotate - bins
params$metannotate_data_filename <- here::here("input_data", "metannotate_annotations_genome_bins.tsv")
params$hmm_info_filename <- here::here("input_data", "hmm_info.tsv")
params$genome_naming_info_filename <- here::here("input_data", "genome_naming_info.tsv")
params$evalue_cutoff <- 1e-40
params$rank_to_summarize <- "family"
params$taxon_to_keep <- "Chlorobia"
params$taxon_keep_rank <- "class"

params$output_filepath <- here::here("plot", "Figure_03_raw.pdf")
##########################

# Function to choose an optimal colour scale to suit the number of entries being used
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

# Function: Helper function to quickly get the number of unique entries in a column of a tibble
count_unique_ranks <- function(stats_table, grouping_var) {
  number_of_ranks <- length(unique(dplyr::pull(ungroup(stats_table), grouping_var)))
  
  return(number_of_ranks)
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
parse_GTDB_taxonomy_entry <- function(taxonomy_entry, resolve = FALSE) {
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

# Function to parse GTDB taxonomy from a vector of taxa into a data frame
parse_GTDB_taxonomy_vector <- function(taxonomy_vector, resolve = FALSE) {
  
  taxonomy <- tibble::as_tibble(as.data.frame(matrix(unlist(
    lapply(taxonomy_vector, function(x) { parse_GTDB_taxonomy_entry(x, resolve = resolve) })), 
    nrow = length(taxonomy_vector), byrow = TRUE), stringsAsFactors = FALSE))
  
  colnames(taxonomy) <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  return(taxonomy)
}

# Function to load GTDB data tables
load_GTDB_data <- function(gtdb_filepath_bac120, gtdb_filepath_ar122, 
                            resolve_ambiguous_taxonomy = FALSE) {
  
  # Load bacteria
  GTDB_data_bacteria <- tibble::as_tibble(read.table(gtdb_filepath_bac120, sep = "\t", header = TRUE,
                                                     stringsAsFactors = FALSE)) %>%
                        dplyr::select(user_genome, classification)
  
  # Load archaea
  # TODO - skip if the archaea input is NA
  GTDB_data_archaea <- tibble::as_tibble(read.table(gtdb_filepath_ar122, sep = "\t", header = TRUE,
                                                     stringsAsFactors = FALSE)) %>%
                       dplyr::select(user_genome, classification)
  # (Note the usage of dplyr::select - must simplify due to some issues with recognizing RED values as numeric vs character)
  
  # Join and simplify
  GTDB_data <- dplyr::bind_rows(GTDB_data_bacteria, GTDB_data_archaea) %>%
               dplyr::rename(genome = user_genome)
  
  # Parse taxonomy
  GTDB_data <- parse_GTDB_taxonomy_vector(GTDB_data$classification, resolve = resolve_ambiguous_taxonomy) %>%
               tibble::add_column(genome = GTDB_data$genome, .before = 1)
  
  return(GTDB_data)
  
}

# Function to load read mapping and assembly stats. Calculates some simple, handy stats.
load_read_stats <- function(mapping_stats_filename, assembled_read_stats_filename,
                            read_recruitment = "assembled_reads") {
  
  # Load read mapping stats
  read_mapping_stats <- tibble::as_tibble(read.table(mapping_stats_filename, 
                                                     sep = "\t", header = TRUE, stringsAsFactors = FALSE))
  
  # Load assembled read stats and add to main table
  assembled_read_stats <- tibble::as_tibble(read.table(assembled_read_stats_filename, 
                                                       sep = "\t", header = TRUE, stringsAsFactors = FALSE))
  read_mapping_stats <- dplyr::left_join(read_mapping_stats, assembled_read_stats, by = "metagenome")
  
  # Calc stats
  read_mapping_stats$percent_assembled_reads <- read_mapping_stats$mapped_assembled_reads /
    read_mapping_stats$metagenome_total_reads * 100
  # Calculate the read recruitment to each genome either relative to total or assembled reads
  if (read_recruitment == "assembled_reads") {
    read_mapping_stats$percent_read_recruitment <- read_mapping_stats$mapped_reads / 
      read_mapping_stats$mapped_assembled_reads * 100
  } else if (read_recruitment == "all_reads") {
    read_mapping_stats$percent_read_recruitment <- read_mapping_stats$mapped_reads / 
      read_mapping_stats$metagenome_total_reads * 100
  } else {
    flog.error(glue::glue("read_recruitment must be either 'assembled_reads' or 'all_reads'. You provided ", 
                          read_recruitment, ". Exiting..."))
    quit(save = "no", status = 1)
  }
  
  # Simplify
  read_mapping_stats <- dplyr::select(read_mapping_stats, genome, metagenome, 
                                      percent_read_recruitment, percent_assembled_reads,
                                      mapped_assembled_reads, metagenome_total_reads)
  
  return(read_mapping_stats)
  
}

# Function to load MetAnnotate data and apply an e-value cutoff
load_metannotate_data <- function(metannotate_data_filename, evalue_cutoff = 1e-40, 
                                  hmm_info_filename) {
  
  # Load data and simplify
  metannotate_data <- tibble::as_tibble(read.table(metannotate_data_filename, sep = "\t", 
                                                   header = TRUE, stringsAsFactors = FALSE)) %>%
    dplyr::rename(genome = Dataset, hmm_name = HMM.Family, evalue = HMM.E.val) %>%
    dplyr::select(genome, hmm_name, evalue, starts_with("Closest"), -Closest.Homolog)
 
  # Fix HMM names and get rid of unwanted HMMs
  hmm_info <- tibble::as_tibble(read.table(hmm_info_filename, sep = "\t", header = TRUE, 
                                           stringsAsFactors = FALSE)) %>%
    dplyr::rename(new_name = HMM.Family) %>%
    dplyr::select(raw_name, new_name)
  metannotate_data <- dplyr::filter(metannotate_data, hmm_name %in% hmm_info$raw_name)
  metannotate_data$hmm_name <- plyr::mapvalues(metannotate_data$hmm_name, from = hmm_info$raw_name,
                                               to = hmm_info$new_name)
  
  # Strip extraneous numerical identifier from genome names
  metannotate_data$genome <- unlist(lapply(metannotate_data$genome, function(x) {
    genome_name <- strsplit(x, split = "_")[[1]]
    genome_name <- glue::glue_collapse(genome_name[1:(length(genome_name)-1)], sep = "_")
    return(genome_name)
  }))
  
  # See the effect of e-value filtration
  counts_unfiltered <- dplyr::group_by(metannotate_data, hmm_name) %>%
                       dplyr::summarise(unfiltered = n())
  
  # Filter by evalue
  metannotate_data <- dplyr::filter(metannotate_data, evalue <= evalue_cutoff)
  
  # See the effect of e-value filtration
  counts_comparison <- dplyr::group_by(metannotate_data, hmm_name) %>%
    dplyr::summarise(filtered = n()) %>%
    dplyr::left_join(counts_unfiltered, by = "hmm_name") %>%
    dplyr::select(hmm_name, unfiltered, filtered)
  flog.info(glue::glue("Filtration effect at e-value cutoff of ", evalue_cutoff, ":"))
  print(counts_comparison)
  
  # Done
  return(metannotate_data)
  
}

# Function to conveniently summarize metannotate output into a single column
# ONLY works if you have exactly two HMMs (fairly strict function for now, just for the Fig. 3 plot)
summarize_metannotate_data <- function(metannotate_data) {
  # Convert the MetAnnotate table into a gene hit summary suitable for the plot
  metannotate_summary <- dplyr::group_by(metannotate_data, genome, hmm_name) %>%
    dplyr::summarise(hits = n())
  metannotate_summary <- reshape2::dcast(metannotate_summary, genome ~ hmm_name, value.var = "hits")
  
  metannotate_summary$gene_summary <- unlist(lapply(1:nrow(metannotate_summary), function(x) {
    gene_1 <- colnames(metannotate_summary[2])
    gene_2 <- colnames(metannotate_summary[3])
    
    if (is.na(metannotate_summary[x,2]) == FALSE && 
        is.na(metannotate_summary[x,3]) == FALSE) {
      return(paste(gene_1, " & ", gene_2, sep = ""))
    } else if (is.na(metannotate_summary[x,2]) == FALSE && 
               is.na(metannotate_summary[x,3]) == TRUE) {
      return(gene_1)
    } else if (is.na(metannotate_summary[x,2]) == TRUE && 
               is.na(metannotate_summary[x,3]) == FALSE) {
      return(gene_2)
    } else if (is.na(metannotate_summary[x,2]) == TRUE && 
               is.na(metannotate_summary[x,3]) == TRUE) {
      return("Not applicable")
    } 
  }))
  
  metannotate_summary <- dplyr::select(metannotate_summary, genome, gene_summary)
  
  return(metannotate_summary)
}

# Function to merge all data
merge_all_data <- function(GTDB_data, read_mapping_stats, metannotate_summary,
         checkM_stats_filepath) {
  
  # Load checkM data
  checkM_stats <- tibble::as_tibble(read.table(checkM_stats_filepath, sep = "\t", comment.char = "",
                                               header = TRUE, stringsAsFactors = FALSE)) %>%
    dplyr::select(genome, completeness, contamination, strain_heterogeneity)
  
  # Fix genome names
  checkM_stats$genome <- tools::file_path_sans_ext(checkM_stats$genome)
  checkM_stats$genome <- gsub(pattern = "[[:punct:]]", x = checkM_stats$genome, replacement = "_")
  GTDB_data$genome <- gsub(pattern = "[[:punct:]]", x = GTDB_data$genome, replacement = "_")
  read_mapping_stats$genome <- gsub(pattern = "[[:punct:]]", x = read_mapping_stats$genome, replacement = "_")
  
  # Join tables
  combined_data <- dplyr::left_join(read_mapping_stats, GTDB_data, by = "genome") %>%
    dplyr::left_join(checkM_stats, by = "genome") %>%
    dplyr::left_join(metannotate_summary, by = "genome")
  
  # TODO
  # write.table(combined_data, file = "test.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  
  return(combined_data)
  
}

# Function to make data frame for panel C (overall taxonomy)
panel_C_prep <- function(combined_data, cutoff = 0.01, rank_to_summarize = "family") {
  # HARD-CODED
  tax_ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  
  # Filter out low coverage genomes
  cutoff <- cutoff * 100
  combined_data_filtered <- dplyr::filter(combined_data, percent_read_recruitment >= cutoff)
  flog.info(glue::glue("Filtering to ", cutoff, "% or higher abundance bins. Reduced from ", 
                       length(unique(dplyr::pull(combined_data, genome))), " to ", 
                       length(unique(dplyr::pull(combined_data_filtered, genome))), " bins."))
  
  # Collapse to taxonomy
  group_by_vars <- c("metagenome", tax_ranks[1:match(x = rank_to_summarize, table = tax_ranks)])
  combined_data_collapsed <- dplyr::group_by_at(combined_data_filtered, group_by_vars) %>%
          dplyr::summarise(percent_read_recruitment = sum(percent_read_recruitment, na.rm = TRUE))
  
  # Calculate stats
  pre_collapse <- count_unique_ranks(combined_data_filtered, "genome")
  post_collapse <- count_unique_ranks(combined_data_collapsed, rank_to_summarize)
  flog.info(glue::glue("Summarized table to rank '", rank_to_summarize, "'. Collapsed from ", 
                       pre_collapse, " genome bins to ", post_collapse, " taxonomic groups"))
  
  if (rank_to_summarize == "genus") {
    # Add family to genus name to make less ambiguous
    combined_data_collapsed$plotting_name <- paste(combined_data_collapsed$genus, " (f__", 
                                                   combined_data_collapsed$family, ")", sep = "")
  } 
  else {
    combined_data_collapsed$plotting_name <- dplyr::pull(combined_data_collapsed, rank_to_summarize)

    # Quick and dirty method to clarify **most** "Unresolved" related issues
    combined_data_collapsed$plotting_name <- unlist(lapply(1:nrow(combined_data_collapsed), function(x) {
      plot_name <- dplyr::pull(combined_data_collapsed, plotting_name)[x]
      if (plot_name == "Unresolved") {
        # Find the next rank up
        rank_colnum <- match(x = rank_to_summarize, table = colnames(combined_data_collapsed))
        prev_rank <- dplyr::pull(combined_data_collapsed, rank_colnum - 1)[x]
        new_name <- paste(plot_name, "_", prev_rank, sep = "")
        if (new_name == "Unresolved_Unresolved") {
          flog.warn(glue::glue("At least one of the plotting names is 'Unresolved_Unresolved'. If you have more
                               than one entry like this, your plot will not be accurate. Check the output table
                               to the user to see if this is the case."))
        } else {
          flog.info(glue::glue("Replaced ambiguous name ", plot_name, " with ", new_name))
        }
        return(new_name)
      } else {
        return(plot_name)
      }
      }))

    }
  
  combined_data_collapsed$plot_panel <- "C"
  
  return(combined_data_collapsed)
}

### Plot for specific genes
panel_A_B_prep <- function(combined_data, cutoff = 0.001, 
                           taxon_to_keep = "Chlorobia", taxon_keep_rank = "class") {
  
  # Reduce to genomes that had a hit to one of the genes or is the search taxon
  combined_data_searched <- dplyr::filter(combined_data, is.na(gene_summary) == FALSE | 
                                   (!!rlang::sym(taxon_keep_rank)) == taxon_to_keep)
  
  # Filter out low coverage genomes
  cutoff <- cutoff * 100
  combined_data_filtered <- dplyr::filter(combined_data_searched, percent_read_recruitment >= cutoff)
  flog.info(glue::glue("Filtering to ", cutoff, "% or higher abundance bins. Reduced from ",
                       length(unique(dplyr::pull(combined_data_searched, genome))), " to ",
                       length(unique(dplyr::pull(combined_data_filtered, genome))), " bins."))
  
  # Add more info to bin names
  checkm_summary <- paste(combined_data_filtered$completeness, "/", combined_data_filtered$contamination, "/",
                          combined_data_filtered$strain_heterogeneity, sep = "")
  taxonomy_summary <- paste("f__", combined_data_filtered$family, 
                            "; g__", combined_data_filtered$genus, sep = "")
  combined_data_filtered$plotting_name <- paste(combined_data_filtered$genome, " (", checkm_summary, 
                                         "); ", taxonomy_summary, sep = "")
  
  # Give a special name to organisms in the taxon_to_keep
  combined_data_filtered$plot_panel <- unlist(lapply(dplyr::pull(combined_data_filtered, taxon_keep_rank), 
                                                     function(x) {
    if (x == taxon_to_keep) {
      return("A")
    } else {
      return("B")
    }
  }))
  
  # Change NA to "Not applicable" before merging with panel C
  combined_data_filtered$gene_summary <- unlist(lapply(combined_data_filtered$gene_summary, function(x) {
    if (is.na(x) == TRUE) {
      return("Not applicable")
    } else {
      return(x)
    }
  }))
  
  return(combined_data_filtered)
  
}

main <- function(params) {
  # Function to load GTDB data tables
  GTDB_data <- load_GTDB_data(gtdb_filepath_bac120 = params$gtdb_filepath_bac120, 
                              gtdb_filepath_ar122 = params$gtdb_filepath_ar122, 
                              resolve_ambiguous_taxonomy = FALSE)
    
  read_mapping_stats <- load_read_stats(mapping_stats_filename = params$mapping_stats_filename,
                                        assembled_read_stats_filename = params$assembled_read_stats_filename,
                                        read_recruitment = "assembled_reads")
      
  metannotate_summary <- load_metannotate_data(metannotate_data_filename = params$metannotate_data_filename, 
                                            evalue_cutoff = params$evalue_cutoff, 
                                            hmm_info_filename = params$hmm_info_filename) %>%
                      summarize_metannotate_data()
  
  combined_data <- merge_all_data(GTDB_data = GTDB_data, read_mapping_stats = read_mapping_stats, 
                                  metannotate_summary = metannotate_summary,
                                  checkM_stats_filepath = params$checkM_stats_filepath)
 
  # Make genome and metagenome names plotting ready; also add order to metagenomes
  genome_naming <- tibble::as_tibble(read.table(params$genome_naming_info_filename, sep = "\t",
                                                header = TRUE, stringsAsFactors = FALSE))
  combined_data$genome <- plyr::mapvalues(combined_data$genome, from = genome_naming$raw_name,
                                          to = genome_naming$new_name)
  metagenome_naming <- tibble::as_tibble(read.table(params$metagenome_naming_info_filename, sep = "\t",
                                                    header = TRUE, stringsAsFactors = FALSE))
  combined_data$metagenome <- plyr::mapvalues(combined_data$metagenome, from = metagenome_naming$raw_name,
                                          to = metagenome_naming$new_name) %>%
                              factor(levels = metagenome_naming$new_name, ordered = TRUE)
  
  # Panels A and B
  panel_A_B_data <- panel_A_B_prep(combined_data, cutoff = 0.001, 
                             taxon_to_keep = "Chlorobia", taxon_keep_rank = "class")
  # Panel C
  panel_C_data <- panel_C_prep(combined_data, cutoff = 0.01, rank_to_summarize = "family") %>%
    tibble::add_column(gene_summary = "Not applicable")
  
  
  # Combine
  plotting_table <- dplyr::bind_rows(panel_A_B_data, panel_C_data)
  
  # MANUALLY assign some factor order
  plotting_table$gene_summary <- factor(plotting_table$gene_summary, levels =
                                              c("cyc2 & dsrA", "cyc2", "dsrA", "Not applicable"),
                                              ordered = TRUE)
  
  
  # Add lake info
  # Add more info for sorting in the plot
  plotting_table$Lake <- unlist(lapply(plotting_table$genome, function(x) {
    if (grepl(pattern = "^L227", x = x)) {
      return("L227")
    } else if (grepl(pattern = "^L442", x = x)) {
      return("L442")
    } else if (grepl(pattern = "^L304", x = x)) {
      return("L304")
    } else {
      return(NA)
    }
  }))
  
  # Sort meaningfully
  plotting_table <- dplyr::arrange(plotting_table, gene_summary, domain, phylum, class, order, family, genus, 
                                        Lake, genome, plotting_name) 
  plotting_table$plotting_name <- factor(plotting_table$plotting_name, 
                                               levels = rev(unique(plotting_table$plotting_name)),
                                               ordered = TRUE)
  
  # Assign colours
  colour_scale <- choose_discrete_colour_scale(length(unique(plotting_table$gene_summary)))
  # Make the last entry grey
  colour_scale[length(colour_scale)] <- "#a6a6a6"
  
  # Print off tables for the user to be aware of
  plotting_table_print <- unique(dplyr::select(plotting_table, domain, phylum, class,
                                              order, family, genus, genome, plot_panel, gene_summary))
  write.table(plotting_table_print, file = "plotting_name_arrangement_raw.tsv", sep = "\t", 
              row.names = FALSE, col.names = TRUE)

  plot_y_axis <- "percent_read_recruitment"
  plotting_table$labels <- round(dplyr::pull(plotting_table, plot_y_axis), digits = 1)
  
  combined_bins_plot <- ggplot(plotting_table, aes_string(x = "metagenome", y = "plotting_name")) +
    # geom_bar(aes_string(weight = plot_y_axis, fill = taxonomic_rank)) +
    geom_point(aes_string(size = plot_y_axis, fill = "gene_summary"), 
               shape = 21, alpha = 0.6) +
    geom_text(aes(label = labels), size = 2) +
    facet_grid(plot_panel ~ ., scales = "free", space = "free") +
    theme_bw() +
    theme(axis.title = element_text(size = 16), axis.text = element_text(size = 10, colour = "black"),
          strip.text = element_blank(), strip.background = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
          # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.ticks = element_line(size = 0.5, colour = "black"), axis.line = element_line(size = 0.5),
          legend.text = element_text(size = 12), legend.title = element_text(size = 12, face = "bold"),
          legend.key = element_rect(colour = "transparent"), legend.key.size = unit(6, "mm"),
          legend.spacing = unit(1, "mm"), legend.box.just = "left",
          panel.grid = element_line(colour = "#999999", size = 0.1)) +
    scale_fill_manual(values = colour_scale) +
    guides(fill = guide_legend(title = "Genes present in bin", override.aes = list(size = 10)),
           size = FALSE) +
           # size = guide_legend(title = "Mapped assembled reads (%)")) +
    xlab("Sample") +
    scale_size_continuous(range = c(3,40), breaks = c(5, 10, 20, 40, 60)) +
    ylab("Genome bin")
  
  # Bar plot of assembly stats
  assembly_stats_table <- unique(dplyr::select(combined_data, metagenome, 
                                               mapped_assembled_reads, metagenome_total_reads))
  assembly_stats_table$unassembled_reads <- assembly_stats_table$metagenome_total_reads - 
    assembly_stats_table$mapped_assembled_reads
  assembly_stats_table <- dplyr::select(assembly_stats_table, -metagenome_total_reads) %>%
                          dplyr::rename(mapped_reads = mapped_assembled_reads, 
                                        unmapped_reads = unassembled_reads)
  assembly_stats_table <- tibble::as_tibble(reshape2::melt(assembly_stats_table, id.vars = "metagenome", 
                                                           variable.name = "count_type", 
                                                           value.name = "read_count"))
  assembly_stats_table_factors <- tibble::tibble(raw_name = c("mapped_reads", "unmapped_reads"),
                                                 plotting_name = c("Assembled", "Unassembled"))
  assembly_stats_table$count_type <- plyr::mapvalues(as.character(assembly_stats_table$count_type),
                                                     from = assembly_stats_table_factors$raw_name,
                                                     to = assembly_stats_table_factors$plotting_name)
  assembly_stats_table$count_type <- factor(assembly_stats_table$count_type, 
                                            levels = rev(assembly_stats_table_factors$plotting_name), ordered = TRUE)
  assembly_stats_plot <- ggplot(assembly_stats_table, aes(x = metagenome)) +
    geom_bar(aes(weight = read_count, fill = count_type)) +
    scale_fill_manual(values = c("#a6a6a6", "#2857e4")) +
    scale_y_continuous(labels = scales::comma) +
    guides(fill = guide_legend(title = "Metagenomic reads")) +
    theme_bw() +
    theme(axis.title = element_text(size = 16), axis.text = element_text(size = 10, colour = "black"),
          strip.text = element_text(size = 11, face = "italic"), strip.background = element_rect(fill = "#e6e6e6"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.ticks = element_line(size = 0.5, colour = "black"), axis.line = element_line(size = 0.5),
          legend.text = element_text(size = 12), legend.title = element_text(size = 12, face = "bold"),
          legend.key = element_rect(colour = "transparent"), legend.key.size = unit(6, "mm"),
          legend.spacing = unit(1, "mm"), legend.box.just = "left",
          panel.grid = element_blank()) +
    xlab("Sample") +
    ylab("Read count")
  
  pdf(file = params$output_filepath, width = 320/25.4, height = 405/25.4, onefile = FALSE)
  # For 'onefile', see https://stackoverflow.com/a/12483347 (accessed 181210)
  egg::ggarrange(plots = list(combined_bins_plot, assembly_stats_plot), ncol = 1, heights = c(12,1))
  dev.off()
  
}

main(params)
