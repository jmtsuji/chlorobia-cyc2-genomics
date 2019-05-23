# Figure_03_plotter.R
# Generates Figure 3 for the Chlorobi cyc2 paper -- environmental abundance of MAGs
# Copyright Jackson M. Tsuji, 2018 (Neufeld Research Group)
# See required packages in the library header here as well as in the headers of the two dependency scripts below

##### Load libraries and dependency scripts
library(here)
# Some functions are used from a script I initially developed for a different purpose for internal lab use
source(here::here("Figure_03_library.R"))
# Some functions from another repo are also used for parsing MetAnnotate info
source("https://github.com/jmtsuji/metannotate-analysis/releases/download/v0.9.5/metannotate_barplots.R")
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
params$dataset_info_read_mapping_filename <- here::here("input_data", "dataset_info_read_mapping.tsv")
params$assembled_read_stats_filename <- here::here("input_data", "assembled_read_stats.tsv")
params$checkM_stats_filepath <- here::here("input_data", "checkm_stats_reduced.tsv")
params$clarify_unresolved_taxa <- FALSE

# MetAnnotate - bins
params$metannotate_bin_filename <- here::here("input_data", "metannotate_annotations_genome_bins.tsv.gz")
params$hmm_info_bins_filename <- here::here("input_data", "hmm_info.tsv")
params$dataset_info_bins_filename <- here::here("input_data", "dataset_info_metannotate.tsv")
params$evalue_cutoff_bins <- 1e-40
params$rank_to_summarize <- "family"
params$taxon_to_keep <- "Chlorobia"
params$taxon_keep_rank <- "class"

params$output_filepath <- here::here("Figure_03_raw.pdf")
##########################


### Load the data - GTDB and read mapping
load_and_process_GTDB_read_mapping_data <- function(mapping_stats_filename, gtdb_filepath_bac120, gtdb_filepath_ar122,
                                                    genome_info_filename, checkM_stats_filepath,
                                                    metagenome_info_filename, assembled_read_stats_filename,
                                                    drep_summary_filename, write_table = FALSE) {
  GTDB_stats_data <- tibble::as.tibble(read_GTDB_abundance_data(input_bin_mapping_stats_filename = mapping_stats_filename,
                                                                input_gtdb_classification_bac120_filename = gtdb_filepath_bac120, 
                                                                input_gtdb_classification_ar122_filename = gtdb_filepath_ar122,
                                                                clarify_unresolved_taxa = FALSE))
  
  # Add on checkM stats
  checkM_stats <- tibble::as.tibble(read.table(checkM_stats_filepath, 
                                             sep = "\t", header = TRUE, stringsAsFactors = FALSE))
  checkM_stats <- dplyr::select(checkM_stats, genome, completeness, contamination, strain_heterogeneity)
  checkM_stats$genome <- tools::file_path_sans_ext(checkM_stats$genome)
  
  GTDB_stats_data <- dplyr::left_join(GTDB_stats_data, checkM_stats, by = "genome")
  
  # Clean up GTDB genome names to match MetAnnotate
  genome_naming <- tibble::as.tibble(read.table(genome_info_filename, 
                                                sep = "\t", header = TRUE, stringsAsFactors = FALSE))
  genome_naming$gtdb <- unlist(lapply(genome_naming$raw_name, function(x) {
    y <- strsplit(x, split = "_")[[1]]
    y <- y[1:(length(y)-1)]
    y <- glue::glue_collapse(y, sep = "_")
    return(y)
  }))
  GTDB_stats_data$genome <- gsub(pattern = "[[:punct:]]", replacement = "_", x = GTDB_stats_data$genome)
  GTDB_stats_data <- dplyr::filter(GTDB_stats_data, genome %in% genome_naming$gtdb)
  GTDB_stats_data$genome <- plyr::mapvalues(x = GTDB_stats_data$genome, from = genome_naming$gtdb, to = genome_naming$Dataset)
  
  # # Order genomes by taxonomy
  # GTDB_stats_data <- dplyr::arrange(GTDB_stats_data, domain, phylum, class, order, family, genus, species, genome)
  # GTDB_stats_data$genome <- factor(GTDB_stats_data$genome, levels = unique(GTDB_stats_data$genome), ordered = TRUE)
  # 
  # Add on assembled read stats; calculate percent read recruitment for assembled reads
  assembly_stats <- tibble::as_tibble(read.table(assembled_read_stats_filename, 
                                                 sep = "\t", header = TRUE, stringsAsFactors = FALSE))
  GTDB_stats_data <- dplyr::left_join(GTDB_stats_data, assembly_stats, by = "metagenome")
  GTDB_stats_data$percent_read_recruitment_assembled <- GTDB_stats_data$mapped_reads /
    GTDB_stats_data$mapped_assembled_reads * 100
  GTDB_stats_data$percent_assembled_reads <- GTDB_stats_data$mapped_assembled_reads /
    GTDB_stats_data$metagenome_total_reads * 100
  
  # Change dataset names to match MetAnnotate
  metagenome_info <- tibble::as.tibble(read.table(metagenome_info_filename, sep = "\t", header = TRUE, 
                                               stringsAsFactors = FALSE, comment.char = "", quote = ""))
  GTDB_stats_data <- dplyr::filter(GTDB_stats_data, (metagenome %in% metagenome_info$raw_name))
  GTDB_stats_data$metagenome <- plyr::mapvalues(x = GTDB_stats_data$metagenome, from = metagenome_info$raw_name, 
                                                to = metagenome_info$Dataset)
  GTDB_stats_data$metagenome <- factor(GTDB_stats_data$metagenome, levels = unique(metagenome_info$Dataset), 
                                       ordered = TRUE)
  
  if (write_table == TRUE) {
    write.table(file = "bin_colouring_info.tsv", x = unique(dplyr::select(GTDB_stats_data, 
      genome, domain, phylum, class, order, family, genus, species, completeness, contamination, strain_heterogeneity)), 
      sep = "\t", row.names = FALSE, col.names = TRUE)
  }
    
  return(GTDB_stats_data)
}


### Attach MetAnnotate bin data onto the GTDB/read mapping table
attach_read_mapping_and_gene_presence_absense_info <- function(metannotate_bin_filename, hmm_info_filename, dataset_info_filename, 
                                                               evalue_cutoff, GTDB_stats_data) {
  metannotate_bin_data <- read_metannotate_tibble(metannotate_bin_filename)
  
  ### Load the data - MetAnnotate, bins
  # # ASSUMES you have already run this
  # setup_templates <- create_setup_templates(metannotate_bin_data, write_tables = TRUE)
  # flog.info(glue::glue("Wrote setup templates to 'hmm_info_template.tsv' and 'dataset_info_template.tsv. "))
  
  # Map on values
  flog.info("Loading user-provided HMM and dataset naming information")
  metannotate_bin_data <- map_naming_information(metannotate_bin_data, hmm_info_filename, dataset_info_filename)
  
  # But remove sorting to keep compatible with GTDB, for the Dataset column
  metannotate_bin_data$Dataset <- as.character(metannotate_bin_data$Dataset)
  
  # Filter to e-value cutoff
  # hits_original <- tibble::as.tibble(summarize_total_reads_all_genes(metannotate_bin_data))
  metannotate_bin_data <- dplyr::filter(metannotate_bin_data, HMM.E.val <= evalue_cutoff)
  
  # This table can ultimately be used for presence/absense of genes
  gene_presence_absence_wide <- tibble::as.tibble(summarize_total_reads_all_genes(metannotate_bin_data))
  gene_presence_absence_long <- tibble::as.tibble(summarize_total_reads_all_genes(metannotate_bin_data, format = "long"))
  
  # Bind on the gene presence/absence info
  gene_presence_absence_long <- dplyr::rename(gene_presence_absence_long, genome = Dataset)
  gene_presence_absence_wide <- dplyr::rename(gene_presence_absence_wide, genome = Dataset)
  GTDB_stats_data <- dplyr::left_join(GTDB_stats_data, gene_presence_absence_wide, by = "genome")
  
  output_list <- list(GTDB_stats_data, gene_presence_absence_long)
  names(output_list) <- c("GTDB_stats_data", "gene_counts_long")
  return(output_list)
  
}


### Overall taxonomy plot
make_overall_plotting_table <- function(GTDB_stats_data_list, cutoff = 0.01, rank_to_summarize = "family") {
  GTDB_stats_data <- GTDB_stats_data_list$GTDB_stats_data
  
  # Filter out low coverage genomes
  cutoff <- cutoff * 100
  GTDB_stats_data_filtered <- dplyr::filter(GTDB_stats_data, percent_read_recruitment_assembled >= cutoff)
  flog.info(glue::glue("Filtering to ", cutoff, "% or higher abundance bins. Reduced from ", 
                       length(unique(dplyr::pull(GTDB_stats_data, genome))), " to ", 
                       length(unique(dplyr::pull(GTDB_stats_data_filtered, genome))), " bins."))
  
  # HARD-CODED
  tax_ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  
  # Sum the stats -- unfortunately, standard deviation is lost here
  group_by_vars <- c("metagenome", tax_ranks[1:match(x = rank_to_summarize, table = tax_ranks)])
  stats_table_grp <- dplyr::group_by_at(GTDB_stats_data_filtered, group_by_vars)
  stats_table_summ <- dplyr::summarise(stats_table_grp, coverage_mean = sum(coverage_mean, na.rm = TRUE),
                                       rpkb = sum(rpkb, na.rm = TRUE), percent_read_recruitment = 
                                         sum(percent_read_recruitment, na.rm = TRUE), percent_read_recruitment_assembled = 
                                         sum(percent_read_recruitment_assembled, na.rm = TRUE))
  
  # Calculate stats
  pre_collapse <- count_unique_ranks(GTDB_stats_data_filtered, "genome")
  post_collapse <- count_unique_ranks(stats_table_summ, rank_to_summarize)
  flog.info(glue::glue("Summarized table to rank '", rank_to_summarize, "'. Collapsed from ", 
                       pre_collapse, " genome bins to ", post_collapse, " taxonomic groups"))
  
  if (rank_to_summarize == "genus") {
    # Add family to genus name to make less ambiguous
    stats_table_summ$plotting_name <- paste(stats_table_summ$genus, " (f__", stats_table_summ$family, ")", sep = "")
  } else {
    stats_table_summ$plotting_name <- dplyr::pull(stats_table_summ, rank_to_summarize)
    
    # Quick and dirty method to clarify **most** "Unresolved" related issues
    stats_table_summ$plotting_name <- unlist(lapply(1:nrow(stats_table_summ), function(x) {
      plot_name <- dplyr::pull(stats_table_summ, plotting_name)[x]
      if (plot_name == "Unresolved") {
        # Find the next rank up
        rank_colnum <- match(x = rank_to_summarize, table = colnames(stats_table_summ))
        prev_rank <- dplyr::pull(stats_table_summ, rank_colnum - 1)[x]
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
  
  # Save a simplified version to map more info onto later
  stats_table_overall_plot <- dplyr::select(dplyr::ungroup(stats_table_summ), metagenome,
                                            tax_ranks[1:match(x = rank_to_summarize, table = tax_ranks)], 
                                            plotting_name, percent_read_recruitment_assembled)
  stats_table_overall_plot$type <- "C_overall"
  
  return(stats_table_overall_plot)
}

### Plot for specific genes
make_gene_plotting_table <- function(GTDB_stats_data_list, cutoff = 0.001,
                                     taxon_to_keep = "Chlorobia", taxon_keep_rank = "class") {
  GTDB_stats_data <- GTDB_stats_data_list$GTDB_stats_data
  
  # Filter out low coverage genomes
  cutoff <- cutoff * 100
  stats_table_genes <- dplyr::filter(GTDB_stats_data, percent_read_recruitment_assembled >= cutoff)
  flog.info(glue::glue("Filtering to ", cutoff, "% or higher abundance bins. Reduced from ",
                       length(unique(dplyr::pull(GTDB_stats_data, genome))), " to ",
                       length(unique(dplyr::pull(stats_table_genes, genome))), " bins."))
  
  
  # Add more info to bin names
  checkm_summary <- paste(stats_table_genes$completeness, "/", stats_table_genes$contamination, "/",
                          stats_table_genes$strain_heterogeneity, sep = "")
  taxonomy_summary <- paste("f__", stats_table_genes$family, 
                            "; g__", stats_table_genes$genus, sep = "")
  stats_table_genes$genome_long <- paste(stats_table_genes$genome, " (", checkm_summary, 
                                         "); ", taxonomy_summary, sep = "")
  
  # Generate gene summary
  genes_to_plot <- as.character(unique(dplyr::pull(GTDB_stats_data_list$gene_counts_long, gene)))
  gene_info <- dplyr::select(stats_table_genes, genes_to_plot)
  
  stats_table_genes$gene_summary <- unlist(lapply(1:nrow(gene_info), function(x) {
    gene_row <- gene_info[x,]
    gene_summ <- glue::glue_collapse(unlist(lapply(1:ncol(gene_row), function(y) {
      if (is.na(dplyr::pull(gene_row, y)) == TRUE) {
        flog.debug(glue::glue(x, ": no ", y))
      } else {
        return(colnames(gene_row)[y])
      }
    })), sep = " & ")
    
    if (length(gene_summ) == 0) {
      gene_summ <- "Not applicable"
    }
    
    return(gene_summ)
  }))
  
  # Keep the requested taxon, BUT omit any other rows lacking all the genes
  # for sym, see https://stackoverflow.com/a/24569597 subcomment (accessed 181211)
  stats_table_genes <- dplyr::filter(stats_table_genes, gene_summary != "Not applicable" | (!!rlang::sym(taxon_keep_rank)) == taxon_to_keep)
  
  # Make a simplified version for combining later
  stats_table_genes_simple <- dplyr::select(ungroup(stats_table_genes), metagenome, domain, phylum, class, order, 
                                            family, genus, genome_long, percent_read_recruitment_assembled, 
                                            # genes_to_plot)
                                            gene_summary)
  stats_table_genes_simple <- dplyr::rename(stats_table_genes_simple, plotting_name = genome_long)
  
  # Give a special name to organisms in the taxon_to_keep
  stats_table_genes_simple$type <- unlist(lapply(1:nrow(stats_table_genes_simple), function(x) {
    rank_name <- dplyr::pull(stats_table_genes_simple, taxon_keep_rank)[x]
    if (rank_name == taxon_to_keep) {
      return("A_unique_taxon")
    } else {
      return("B_other")
    }
  }))
  
  return(stats_table_genes_simple)
  
}

main <- function(params) {
  # Load GTDB data and correct genome names
  GTDB_stats_data <- load_and_process_GTDB_read_mapping_data(mapping_stats_filename = params$mapping_stats_filename, 
                                                             gtdb_filepath_bac120 = params$gtdb_filepath_bac120, 
                                                             gtdb_filepath_ar122 = params$gtdb_filepath_ar122,
                                                             genome_info_filename = params$dataset_info_bins_filename,
                                                             metagenome_info_filename = params$dataset_info_read_mapping_filename,
                                                             assembled_read_stats_filename = params$assembled_read_stats_filename,
                                                             checkM_stats_filepath = params$checkM_stats_filepath)
  
  # Add on MetAnnotate bin data for presence/absense of genes
  GTDB_stats_data_list <- attach_read_mapping_and_gene_presence_absense_info(
                                                     metannotate_bin_filename = params$metannotate_bin_filename, 
                                                     hmm_info_filename = params$hmm_info_bins_filename, 
                                                     dataset_info_filename = params$dataset_info_bins_filename, 
                                                     evalue_cutoff = params$evalue_cutoff_bins, 
                                                     GTDB_stats_data = GTDB_stats_data)
  GTDB_stats_data <- GTDB_stats_data_list$GTDB_stats_data
  
  # Make the base plotting tables for the overall genus plot and gene-specific plot
  overall_plotting_table <- make_overall_plotting_table(GTDB_stats_data_list, cutoff = 0.01, 
                                                        rank_to_summarize = params$rank_to_summarize)
  gene_plotting_table <- make_gene_plotting_table(GTDB_stats_data_list, cutoff = 0.001, 
                                                  taxon_to_keep = params$taxon_to_keep, 
                                                  taxon_keep_rank = params$taxon_keep_rank)
  
  # Prepare for merging
  overall_plotting_table$gene_summary <- "Not applicable"
  
  ### Make combined plot
  # # Map the gene occurrance patterns of bins to their corresponding genus on the genus plot
  # unique(dplyr::select(gene_plotting_table, family, genus, gene_summary))
  
  stats_table_combined <- dplyr::bind_rows(overall_plotting_table, gene_plotting_table)
  
  # Add more info for sorting in the plot
  stats_table_combined$Lake <- unlist(lapply(1:nrow(stats_table_combined), function(x) {
    if (grepl(pattern = "^L227", x = stats_table_combined$plotting_name[x])) {
      return("L227")
    } else if (grepl(pattern = "^L442", x = stats_table_combined$plotting_name[x])) {
      return("L442")
    } else if (grepl(pattern = "^L304", x = stats_table_combined$plotting_name[x])) {
      return("L304")
    } else {
      return(NA)
    }
  }))
  
  stats_table_combined$gene_summary <- factor(stats_table_combined$gene_summary, levels =
                                                c("cyc2 & dsrA", "cyc2", "dsrA", "Not applicable"),
                                              ordered = TRUE)
  
  # Sort meaningfully
  stats_table_combined <- dplyr::arrange(stats_table_combined, gene_summary, domain, phylum, class, order, family, genus, 
                                            Lake, plotting_name) 
  
  # Now arrange the plotting name based on the above determined sort order
  stats_table_combined$plotting_name <- factor(stats_table_combined$plotting_name, 
                                               levels = rev(unique(stats_table_combined$plotting_name)),
                                               ordered = TRUE)
  
  colour_scale <- choose_discrete_colour_scale(length(unique(stats_table_combined$gene_summary)))
  
  # Make the last entry grey
  colour_scale[length(colour_scale)] <- "#a6a6a6"
  
  # Print off tables for the user to be aware of
  stats_table_combined_print <- unique(dplyr::select(stats_table_combined, domain, phylum, class,
                                              order, family, genus, plotting_name, type, gene_summary))
  write.table(stats_table_combined_print, file = "plotting_name_arrangement_raw.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)

  plot_y_axis <- "percent_read_recruitment_assembled"
  stats_table_combined$labels <- round(dplyr::pull(stats_table_combined, plot_y_axis), digits = 1)
  
  combined_bins_plot <- ggplot(stats_table_combined, aes_string(x = "metagenome", y = "plotting_name")) +
    # geom_bar(aes_string(weight = plot_y_axis, fill = taxonomic_rank)) +
    geom_point(aes_string(size = plot_y_axis, fill = "gene_summary"), 
               shape = 21, alpha = 0.6) +
    geom_text(aes(label = labels), size = 2) +
    facet_grid(type ~ ., scales = "free", space = "free") +
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
  # print(combined_bins_plot)
  # ggsave(file = "comb_plot_vs1.pdf", width = 350, height = 250, units = "mm")
  
  ### Bar plot of assembly stats
  assembly_stats_table <- unique(dplyr::select(GTDB_stats_data, metagenome, metagenome_total_reads,mapped_assembled_reads))
  assembly_stats_table$unassembled_reads <- assembly_stats_table$metagenome_total_reads - assembly_stats_table$mapped_assembled_reads
  assembly_stats_table <- dplyr::select(assembly_stats_table, -metagenome_total_reads)
  assembly_stats_table <- dplyr::rename(assembly_stats_table, mapped_reads = mapped_assembled_reads, 
                                        unmapped_reads = unassembled_reads)
  assembly_stats_table <- tibble::as_tibble(reshape2::melt(assembly_stats_table, id.vars = "metagenome", variable.name = "count_type", 
                                         value.name = "read_count"))
  assembly_stats_table_factors <- tibble::tibble(raw_name = c("mapped_reads", "unmapped_reads"),
                                                 plotting_name = c("Assembled", "Unassembled"))
  assembly_stats_table$count_type <- as.character(assembly_stats_table$count_type)
  assembly_stats_table$count_type <- plyr::mapvalues(assembly_stats_table$count_type,
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
  
  # print(assembly_stats_plot)
  
  pdf(file = params$output_filepath, width = 320/25.4, height = 405/25.4, onefile = FALSE)
  # For 'onefile', see https://stackoverflow.com/a/12483347 (accessed 181210)
  egg::ggarrange(plots = list(combined_bins_plot, assembly_stats_plot), ncol = 1, heights = c(12,1))
  dev.off()
  
}

main(params)
