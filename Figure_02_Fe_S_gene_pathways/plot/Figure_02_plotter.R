# Chlorobia cyc2 manuscript, Figure 02 plotter
# Copyright Jackson M. Tsuji, 2019
# Neufeld Research Group

# Load libraries
library(here)
library(futile.logger)
library(tools)
library(glue, warn.conflicts = FALSE)
library(plyr, warn.conflicts = FALSE)
library(tibble, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(ggplot2, warn.conflicts = FALSE)
library(ggtree, quietly = TRUE, warn.conflicts = FALSE)
library(ape, warn.conflicts = FALSE)
library(maps, warn.conflicts = FALSE)
library(phytools, warn.conflicts = FALSE)
library(reshape2, warn.conflicts = FALSE)
library(gridExtra, warn.conflicts = FALSE)
library(egg, warn.conflicts = FALSE)

####################################
# User variables
gene_naming_info_filename <- here::here("input_data", "gene_naming_info.tsv")
Chlorobia_naming_info_filename <- here::here("input_data", "Chlorobia_naming_info.tsv")
query_subject_filename_separator <- "__to__" # unique character string consistently separating 'queryName' and 'subjectName' in the .csv filenames output by BackBLAST
backblast_header_names <- c("qseqid", "sseqid", "pident", "evalue", "qcovhsp", "bitscore")
backblast_directory_path <- here::here("input_data", "backblast_results")
output_table_filepath <- here::here("plot", "Figure_02_plotting_data.tsv")
output_pdf_filepath <- here::here("plot", "Figure_02_raw.pdf")
phylogenetic_tree_filepath <- here::here("input_data", "chlorobia_riboprotein_phylogeny.treefile")
bootstrap_cutoff <- 50
tree_root <- "Ignavibacterium_album_JCM_16511_outgroup"
####################################

# Function to load table of data, add header IDs, and add query/subject name by parsing filename
load_backblast_results_table <- function(backblast_table_filepath, header_names, query_subject_filename_separator) {
  
  # Get subject and query names
  table_filename_base <- basename(tools::file_path_sans_ext(backblast_table_filepath))
  table_filename_split <- strsplit(table_filename_base, split = query_subject_filename_separator)[[1]]
  if (length(table_filename_split) != 2) {
    flog.fatal(glue::glue(basename(backblast_table_filepath), ": splitting filename by '", query_subject_filename_separator, 
                          "' results in ", length(table_filename_split), " entries (looking for 2). ",
                          "Cannot parse query/subject names. Exiting..."))
    quit(status = 1)
  }
  query_name <- table_filename_split[1]
  subject_name <- table_filename_split[2]
  
  if (file.info(backblast_table_filepath)$size == 0) {
    
    # If the file is empty, warn the user and then return an empty tibble
    flog.warn(glue::glue(basename(backblast_table_filepath), ": file is empty."))
    data_table <- tibble::as_tibble(matrix(c(query_name, subject_name), nrow = 1))
    colnames(data_table) <- c("query_name", "subject_name")
    # Later, dplyr::bind_rows will automatically handle assigning the NA's to the additional columns.
    
  } else {
  
    # Read the table and add columns
    data_table <- tibble::as_tibble(read.table(backblast_table_filepath, header = FALSE, sep = ",", 
                                               stringsAsFactors = FALSE, comment.char = ""))
    colnames(data_table) <- header_names
    
    # Add new columns to table
    data_table <- tibble::add_column(data_table, subject_name = subject_name, .before = 1)
    data_table <- tibble::add_column(data_table, query_name = query_name, .before = 1)
     
  }
  
  return(data_table)
  
}

# Loads all BackBLAST tables, combines, and matches to gene naming table entries
process_backblast_table <- function(backblast_directory_path, gene_naming_info_filename) {

  # Load naming tables
  gene_naming_info <- tibble::as_tibble(read.table(gene_naming_info_filename, sep = "\t", header = TRUE, 
                                                   stringsAsFactors = FALSE, comment.char = "#"))
  
  # Get the backblast table filenames
  backblast_table_filenames <- list.files(path = backblast_directory_path, pattern = "*.csv", full.names = TRUE)
  
  # Combine data from all backblast tables
  backblast_data <- lapply(backblast_table_filenames, function(x) { load_backblast_results_table(
    backblast_table_filepath = x, header_names = backblast_header_names, 
    query_subject_filename_separator = query_subject_filename_separator) } )
  backblast_data <- dplyr::bind_rows(backblast_data)
  
  # Before mapping qseqid's and genes, make sure that the two have perfect correspondance, or else identify and remove the extras
  if (isTRUE(gene_naming_info$gene_accession %in% backblast_data$qseqid) == FALSE) {
    # Then some are missing
    missing_IDs <- gene_naming_info$gene_accession[!(gene_naming_info$gene_accession %in% backblast_data$qseqid)]
    missing_IDs_print <- glue::glue_collapse(missing_IDs, sep = ", ")
    flog.warn(glue::glue("Some accessions had no BLAST hits: ", missing_IDs_print, ". Removing these from the final plot."))
    
    # Remove missing entries
    gene_naming_info <- gene_naming_info[(gene_naming_info$gene_accession %in% backblast_data$qseqid),]
  }
  
  # Map on the proper gene names for clarity (and order)
  gene_naming_info_join <- dplyr::select(gene_naming_info, gene_accession, gene_name)
  backblast_data <- dplyr::inner_join(backblast_data, gene_naming_info_join, by = c("qseqid" = "gene_accession"))
  backblast_data$gene_name <- factor(backblast_data$gene_name, levels = gene_naming_info$gene_name, ordered = TRUE)
  
  return(backblast_data)
}
  
# Loads and plots the phylogenetic tree, re-rooting to the desired outgroup and selecting a bootstrap label cutoff
plot_phylogenetic_tree <- function(phylogenetic_tree_filepath) {
  # Load phylogenetic tree
  flog.info(glue::glue("Loading phylogenetic tree '", phylogenetic_tree_filepath, "'"))
  phylo_tree <- ggtree::read.tree(phylogenetic_tree_filepath)
  
  # Re-root the tree
  tip_label_index <- match(x = tree_root, table = phylo_tree$tip.label)
  phylo_tree <- ggtree::reroot(phylo_tree, node = tip_label_index)
  
  # Plot the tree
  flog.info("Generaring the tree plot")
  tree_plot <- ggtree(phylo_tree, size = 1, colour = "black", ladderize = TRUE,
                      branch.length = 0.1) +
    geom_treescale(x = 0.4, y = 10, linesize = 1, fontsize = 3, offset = 0.2, width = 0.2) +
    geom_tiplab(align = TRUE, linetype = "dotted",
                size = 0, offset = 0.1) +
    geom_text2(aes(subset = (grepl(pattern = "^[0-9]+$", x = label) & !(isTip) & as.numeric(label) > bootstrap_cutoff), 
                   label = as.numeric(label)),
               nudge_x = -0.03, nudge_y = 0.4, size = 2.5) +
    scale_y_discrete(expand = c(0,0.6)) # to manually make it correspond to the heatmap in y-coordinates
  
  return(tree_plot)
}

# Plots the backblast heatmap, aligning columns with the tree
plot_backblast_heatmap <- function(backblast_data, tree_plot, Chlorobia_naming_info_filename) {
  # Load the naming table
  Chlorobia_naming_info <- tibble::as_tibble(read.table(Chlorobia_naming_info_filename, sep = "\t", header = TRUE, 
                                                        stringsAsFactors = FALSE))
  
  # Change the backblast names to match the tree names
  flog.info("Matching backblast table to tree tips")
  backblast_data$subject_name <- plyr::mapvalues(backblast_data$subject_name, from = Chlorobia_naming_info$original_name, 
                                                 to = Chlorobia_naming_info$tree_name)
  
  # Get the order of the tips in the tree plot
  tip_order <- dplyr::filter(tree_plot$data, isTip == TRUE)
  tip_order <- tip_order[order(tip_order$y, decreasing = TRUE),]$label
  
  # Check that the exact same genome names exist between the FastANI data and the tree
  if (all.equal(sort(tip_order), sort(unique(backblast_data$subject_name))) == FALSE) {
    flog.warn("backblast labels and tree tips DO NOT MATCH. Results will be unreliable!")
  }
  
  # Order the genomes in the backblast table to correspond to the tree tip order
  backblast_data$subject_name <- factor(backblast_data$subject_name, levels = rev(tip_order), ordered = TRUE)
  
  # Change taxon names to those to appear in the final plot
  backblast_data$subject_name <- plyr::mapvalues(backblast_data$subject_name, from = Chlorobia_naming_info$tree_name, 
                                                 to = Chlorobia_naming_info$plotting_name)
  
  # If there are multiple blast hits for any subject-query combo, and keep the best hit
  # TODO - warn the user?
  backblast_data <- dplyr::ungroup(backblast_data) %>%
                    dplyr::group_by(subject_name, query_name, qseqid) %>%
                    dplyr::top_n(n = 1, wt = pident) %>%
                    dplyr::ungroup()
  
  # Make the plot
  backblast_heatmap <- ggplot(backblast_data, aes(y = subject_name, x = gene_name)) +
    geom_tile(aes(fill = pident)) +
    guides(fill = guide_colourbar(title = "Amino acid \nidentity (%)")) +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.title = element_text(size = 12),
          strip.text = element_text(size = 7), strip.background = element_rect(fill = "#e6e6e6"),
          panel.border = element_rect(colour = "black", size = 1),
          axis.text = element_text(size = 10, colour = "black", face = "italic"), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "italic"),
          axis.ticks = element_line(size = 0.5), axis.line = element_line(colour = "black", size = 0.5),
          legend.text = element_text(size = 10, colour = "black"),
          legend.title = element_text(size = 10, face = "bold"),
          legend.key = element_blank(), legend.key.size = unit(5, "mm")) +
    xlab("") +
    ylab("")
  
  return(backblast_heatmap)
  
}

main <- function() {
  # Load backblast data
  flog.info("Loading backblast tables and combining")
  backblast_data <- process_backblast_table(backblast_directory_path, gene_naming_info_filename)
  
  # This is a nice checkpoint. Export the table as raw plotting data.
  flog.info("Writing plotting data to file")
  write.table(backblast_data, file = output_table_filepath, 
              sep = "\t", row.names = FALSE, col.names = TRUE)
  
  # Make the plots
  flog.info("Making tree plot")
  tree_plot <- plot_phylogenetic_tree(phylogenetic_tree_filepath)
  flog.info("Making backblast heatmap")
  backblast_heatmap <- plot_backblast_heatmap(backblast_data, tree_plot, Chlorobia_naming_info_filename)
  
  # Put everything together and save
  flog.info(glue::glue("Combining tree and heatmap; saving to '", output_pdf_filepath, "'"))
  pdf(file = output_pdf_filepath, width = 14, height = 4.75, onefile = FALSE)
  egg::ggarrange(tree_plot, backblast_heatmap, nrow = 1,
                 widths = c(1.5,2.5))
  dev.off()
   
}

main()
