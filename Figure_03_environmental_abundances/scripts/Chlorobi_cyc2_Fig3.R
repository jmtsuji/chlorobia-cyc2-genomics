# Chlorobi_cyc2_Fig3.R
# Generates Figure 3 for the Chlorobi cyc2 paper -- environmental abundance of genes in metagenomes
# Copyright Jackson M. Tsuji, 2018 (Neufeld Research Group)

## Load libraries and dependency scripts
setwd("/home/jmtsuji/Research_General/PhD/04b_Metagenome_resequencing_F2015/10_ATLAS_re_analysis/09_env_abundance_of_bins/vs2_dRep_scripted/output/summary/")

source("../../scripts/env_bin_abundance_analysis.R")
source("/home/jmtsuji/Research_General/Bioinformatics/02_git/metannotate-analysis/metannotate_barplots.R")
library(egg)
###

### User-defined variables
params <- list()

# GTDB and read mapping
params$gtdb_filename <- "gtdbtk.summary.tsv"
params$mapping_stats_filename <- "bin_read_mapping_summary.tsv"
params$dataset_info_read_mapping_filename <- "analysis_files/dataset_info_read_mapping_vs1.tsv"
params$assembled_read_stats_filename <- "assembled_read_stats_vs2.tsv"
params$drep_summary_filename <- "checkm_stats_reduced_mod.tsv"
params$clarify_unresolved_taxa <- FALSE

# MetAnnotate - bins
params$metannotate_bin_filename <- "metannotate_dRep_bins.tsv"
params$hmm_info_bins_filename <- "analysis_files/hmm_info_bins_vs6.tsv"
params$dataset_info_bins_filename <- "analysis_files/dataset_info_bins_vs1.tsv"
params$evalue_cutoff_bins <- 1e-40
params$rank_to_summarize <- "family"
params$taxon_to_keep <- "Chlorobia"
params$taxon_keep_rank <- "class"

# MetAnnotate - raw reads
params$metannotate_raw_reads_filename <- "all_annotations_metannotate_raw_reads.tsv"
params$hmm_info_raw_reads_filename <- "analysis_files/hmm_info_raw_reads_vs3.tsv"
params$dataset_info_raw_reads_filename <- "analysis_files/dataset_info_raw_reads_vs1.tsv"
params$evalue_cutoff_raw_reads <- 1e-10
params$taxon_raw_reads <- "Family"
params$raw_read_normalizing_HMM <- "rpoB"

params$output_filepath <- "chlorobi_cyc2_fig3_vs2c.pdf"
#####




### Load the data - GTDB and read mapping
load_and_process_GTDB_read_mapping_data <- function(mapping_stats_filename, gtdb_filename, genome_info_filename,
                                                    metagenome_info_filename, assembled_read_stats_filename,
                                                    drep_summary_filename, write_table = FALSE) {
  GTDB_stats_data <- tibble::as.tibble(read_GTDB_abundance_data(mapping_stats_filename, gtdb_filename, 
                                                                clarify_unresolved_taxa = params$clarify_unresolved_taxa))
  
  # Add on checkM stats
  drep_stats <- tibble::as.tibble(read.table(drep_summary_filename, sep = "\t", header = TRUE, stringsAsFactors = FALSE))
  drep_stats <- dplyr::select(drep_stats, genome, completeness, contamination, strain_heterogeneity)
  drep_stats$genome <- tools::file_path_sans_ext(drep_stats$genome)
  
  GTDB_stats_data <- dplyr::left_join(GTDB_stats_data, drep_stats, by = "genome")
  
  # Clean up GTDB genome names to match MetAnnotate
  genome_naming <- tibble::as.tibble(read.table(genome_info_filename, sep = "\t", header = TRUE, stringsAsFactors = FALSE))
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
  assembly_stats <- read.table(assembled_read_stats_filename, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
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

# Load raw read metannotate data
load_and_process_metannotate_raw_read_data <- function(metannotate_raw_reads_filename, hmm_info_filename,
                                                       dataset_info_filename, evalue_cutoff, taxon,
                                                       normalizing_HMM = "rpoB") {
  
  metannotate_data <- read_metannotate_tibble(metannotate_raw_reads_filename)
  
  ### Load the data - MetAnnotate, bins
  # # ASSUMES you have already run this
  # setup_templates <- create_setup_templates(metannotate_data, write_tables = TRUE)
  # flog.info(glue::glue("Wrote setup templates to 'hmm_info_template.tsv' and 'dataset_info_template.tsv. "))
  
  # Map on values
  flog.info("Loading user-provided HMM and dataset naming information")
  metannotate_data <- map_naming_information(metannotate_data, hmm_info_filename, dataset_info_filename)
  
  # Filter to e-value cutoff
  # hits_original <- tibble::as.tibble(summarize_total_reads_all_genes(metannotate_bin_data))
  flog.info(glue::glue("Filtering by e-value cutoff of ", evalue_cutoff))
  metannotate_data_filtered <- filter_by_evalue(metannotate_data, evalue = evalue_cutoff)
  metannotate_data <- metannotate_data_filtered$metannotate_data
  flog.info(glue::glue("Percent change from e-value filtration:"))
  print(metannotate_data_filtered$read_counts$percent_change)
  
  # Collapse the table to the desired taxonomic rank
  flog.info(glue::glue("Collapsing table to taxonomic rank '", taxon, "'"))
  metannotate_data_collapsed <- collapse_metannotate_table_by_taxon(metannotate_data, taxon = taxon)
  
  # Normalize the data by HMM length
  flog.info("Normalizing data")
  metannotate_data_normalized_list <- normalize_collapsed_metannotate_data(metannotate_data_collapsed, 
                                                                           normalizing_HMM = normalizing_HMM)
  flog.info("Total normalized % abundance of analyzed genes compared to the marker gene:")
  print(metannotate_data_normalized_list$total_normalized_hits)
  
  return(metannotate_data_normalized_list)
  
}

# Contains lots of HARD-CODED variables!!
plot_metannotate_raw_read_data <- function(metannotate_raw_read_data_normalized_list,
                                           colour_template_filename = "colours.tsv") {
  # Plot MetAnnotate data - raw reads - as bubble plot
  # HARD-CODED table relating canonical taxonomy values to the colnames of the metannotate table
  taxonomy_naming <- tibble::tibble(taxonomy = c("domain", "phylum", "class", "order", "family", "genus", "species"),
                                    metannotate_colnames = c("Closest.Homolog.Superkingdom", "Closest.Homolog.Phylum",
                                                             "Closest.Homolog.Class", "Closest.Homolog.Order", "Closest.Homolog.Family", 
                                                             "Closest.Homolog.Genus", "Closest.Homolog.Species"))
  
  # Extract list components
  metannotate_data <- metannotate_raw_read_data_normalized_list$metannotate_data_normalized
  hit_totals_wide <- metannotate_raw_read_data_normalized_list$total_normalized_hits
  hit_totals <- reshape2::melt(hit_totals_wide, id.vars = "Dataset",
                               variable.name = "HMM.Family", value.name = "percent_abundance")
  hit_totals$HMM.Family <- factor(hit_totals$HMM.Family, levels = unique(hit_totals$HMM.Family), ordered = TRUE)
  
  # Detect the taxonomy that the data has been collapsed to
  collapse_taxonomy_metannotate <- tail(taxonomy_naming$metannotate_colnames[
    taxonomy_naming$metannotate_colnames %in% colnames(metannotate_data)], n = 1)
  collapse_taxonomy <- taxonomy_naming$taxonomy[
    match(x = collapse_taxonomy_metannotate, table = taxonomy_naming$metannotate_colnames)]
  flog.debug(glue::glue("Plotting input dataframe has been collapsed to the '", collapse_taxonomy, "' level."))
  
  # Subset to the desired top_x cutoff
  metannotate_data <- subset_normalized_metannotate_data(metannotate_data, top_x = 0.01, percent_mode = "within_sample")
  
  # Make or read in a plotting colour table; or generate auto-colours
  plotting_colours <- process_plotting_colours(metannotate_data, colouring_template_filename= colour_template_filename)
  
  # Make the plotting column into an ordered factor based on the plotting_colours order
  metannotate_data[,collapse_taxonomy_metannotate] <- factor(dplyr::pull(metannotate_data, collapse_taxonomy_metannotate),
                                                             levels = unique(dplyr::pull(plotting_colours, collapse_taxonomy_metannotate)),
                                                             ordered = TRUE)
  
  # Determine normalizing HMM for labelling on the plot
  normalizing_HMM <- colnames(hit_totals_wide)[
    match(100, unlist(lapply(2:ncol(hit_totals_wide), function(x) { mean(hit_totals_wide[,x]) }))) + 1 ]
  hit_totals <- dplyr::filter(hit_totals, HMM.Family != normalizing_HMM)
  
  flog.info("Creating the ggplot")
  
  # Make the bubble plot
  metannotate_plot <- ggplot(metannotate_data) +
    facet_grid(HMM.Family ~ ., scales = "free", space = "free_y") +
    theme_bw() +
    theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10, colour = "black"),
          strip.text = element_text(size = 11, face = "italic"), strip.background = element_rect(fill = "#e6e6e6"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.ticks = element_line(size = 0.5, colour = "black"), axis.line = element_line(size = 0.5),
          legend.text = element_text(size = 10, face = "italic"), legend.title = element_text(size = 8),
          legend.key = element_rect(colour = "transparent"), legend.key.size = unit(6, "mm"),
          legend.spacing = unit(1, "mm"), legend.box.just = "left") +
    xlab("Sample")
  
  metannotate_plot <- metannotate_plot +
    geom_point(aes_string(x = "Dataset", y = collapse_taxonomy_metannotate,
                          size = "percent_abundance", fill = "Closest.Homolog.Class"), shape = 21,
               alpha = 0.8) +
    scale_size_continuous(range = c(1,10)) +
    scale_fill_manual(values = 
                        choose_discrete_colour_scale(length(unique(
                          dplyr::pull(metannotate_data, Closest.Homolog.Class))))) +
    theme(axis.text.y = element_text(size = 5, face = "italic")) +
    guides(size = guide_legend(title = paste("Gene hits relative to \n", normalizing_HMM, 
                                             " (%; normalized)", sep = ""), 
                               override.aes = list(fill = "#4d4d4d")),
           fill = guide_legend(title = "Class")) +
    ylab(paste(strsplit(collapse_taxonomy_metannotate, split = ".", fixed = TRUE)[[1]][3],
               " of closest homologue", sep = ""))
  
  print(metannotate_plot)
  
  return(metannotate_plot)
  
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
  stats_table_overall_plot <- dplyr::select(ungroup(stats_table_summ), metagenome,
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
  GTDB_stats_data <- load_and_process_GTDB_read_mapping_data(params$mapping_stats_filename, params$gtdb_filename,
                                                             genome_info_filename = params$dataset_info_bins_filename,
                                                             metagenome_info_filename = params$dataset_info_read_mapping_filename,
                                                             params$assembled_read_stats_filename,
                                                             params$drep_summary_filename)
  
  # Add on MetAnnotate bin data for presence/absense of genes
  GTDB_stats_data_list <- attach_read_mapping_and_gene_presence_absense_info(
                                                     metannotate_bin_filename = params$metannotate_bin_filename, 
                                                     hmm_info_filename = params$hmm_info_bins_filename, 
                                                     dataset_info_filename = params$dataset_info_bins_filename, 
                                                     evalue_cutoff = params$evalue_cutoff_bins, 
                                                     GTDB_stats_data = GTDB_stats_data)
  GTDB_stats_data <- GTDB_stats_data_list$GTDB_stats_data
  
  # # Load metannotate data
  # metannotate_raw_read_data_normalized_list <- load_and_process_metannotate_raw_read_data(params$metannotate_raw_reads_filename, 
  #                                           params$hmm_info_raw_reads_filename, params$dataset_info_raw_reads_filename,
  #                                           params$evalue_cutoff_raw_reads, params$taxon_raw_reads, params$raw_read_normalizing_HMM) 
  # 
  # # Raw read plot
  # metannotate_raw_read_plot <- plot_metannotate_raw_read_data(metannotate_raw_read_data_normalized_list,
  #                                                             colour_template_filename = "colours_vs2.tsv")
  # 
  
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
  # stats_table_combined$gene_summary <- factor(stats_table_combined$gene_summary, levels =
  #                                               c("cyc2 & dsrA & bchL", "cyc2 & bchL", "dsrA & bchL", 
  #                                                 "cyc2", "dsrA", "bchL", "Not applicable"),
  #                                             ordered = TRUE)
  
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
  # colouring_info <- tibble::tibble(order = unique(dplyr::pull(stats_table_combined, order)), colour = colour_scale)
  # write.table(colouring_info, file = "order_colours_raw.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)
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

# main(params)
