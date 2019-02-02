# plot_cyc2_phylogenetic_tree.R
# Copyright Jackson M. Tsuji, 2018
# Generates one section of the cyc2 figure for the Chlorobi cyc2 paper
# N.B., LOTS of hard-coded tweaks for cyc2; don't use arbitrarily for any tree + MSA

#####################################################
## Load required packages: ##########################
library(futile.logger)
library(plyr)
library(tidyverse, quietly = TRUE, warn.conflicts = FALSE)
library(ggtree, quietly = TRUE, warn.conflicts = FALSE)
library(reshape2, warn.conflicts = FALSE)
library(scales, quietly = TRUE)
library(RColorBrewer)
library(egg)
# library(xtable) # maybe an egg dependency?
library(bio3d)
library(ape)
library(phytools)
#####################################################


create_alignment_data_frame <- function(msa_data) {
  # Convert to data frame
  alignment <- tibble::as.tibble(t(msa_data$ali))
  alignment$position <- 1:nrow(alignment)
  alignment <- as.tibble(reshape2::melt(alignment, id.vars = "position", variable.name = "genome", value.name = "residue"))
  alignment$genome <- as.character(alignment$genome)
  alignment <- dplyr::filter(alignment, residue != "-") # Remove dashes
  
  return(alignment)
}

generate_alignment_summary_stats <- function(alignment, msa_data) {
  # Make summary
  alignment_summ <- tibble(position = 1:max(alignment$position))
  alignment_summ$similarity <- bio3d::conserv(msa_data, method = "similarity", sub.matrix = "blosum62")
  alignment_summ$label <- unlist(lapply(alignment_summ$similarity, function(x) {
    if (x == 1) {
      return("*")
    } else {
      return(NA)
    }
  }))
  # Manually add predicted heme-binding locations
  alignment_summ$heme_binding <- logical(length = nrow(alignment_summ))
  alignment_summ$heme_binding[c(56,59,60)] <- "TRUE"
  
  return(alignment_summ)
}

main <- function(params) {
  
  # Read tree
  flog.info("Reading input phylogenetic tree")
  phylo_tree <- ggtree::read.tree(params$input_phylogenetic_tree_filepath)
  phylo_tree <- phytools::midpoint.root(phylo_tree)
  
  # Plot tree
  flog.info("Plotting tree")
  tree_plot <- ggtree(phylo_tree, size = 1.5, colour = "black", ladderize = TRUE,
                      branch.length = 0.1) +
    # geom_hilight(node = 29, fill = "#1B9E77", alpha = 0.6) + # HARD-CODED to capture Chlorobi
    geom_treescale(x = 0.4, y = 10, linesize = 1, fontsize = 3, offset = 0.2, width = 0.2) +
    geom_tiplab(align = TRUE, linetype = "dotted",
                size = 0, offset = 0.1) +
    geom_text2(aes(subset = (grepl(pattern = "^[0-9]+$", x = label) & !(isTip) & as.numeric(label) > params$bootstrap_cutoff), 
                   label = as.numeric(label)),
               nudge_x = -0.05, nudge_y = 0.4, size = 3) +
    scale_y_discrete(expand = c(0,0.6)) # to manually make it correspond to the heatmap in y-coordinates
    
  # Load MSA data and process
  flog.info("Loading MSA")
  msa_data <- bio3d::read.fasta(file = params$input_msa_filename, rm.dup = TRUE, to.dash = TRUE)
  alignment <- create_alignment_data_frame(msa_data)
  alignment_summ <- generate_alignment_summary_stats(alignment, msa_data)
  
  # Add genome order corresponding to tree tip order
  # TODO - check that the exact same rows exist between the alignment data and the tree tips
  tip_order <- dplyr::filter(tree_plot$data, isTip == TRUE)
  tip_order <- tip_order[order(tip_order$y, decreasing = TRUE),]$label
  alignment$genome <- factor(alignment$genome, levels = rev(tip_order), ordered = TRUE)
  
  # Also change names to more readable names
  plotting_names <- read.table(params$plotting_names_table_filepath, sep = "\t", header = TRUE,
                               stringsAsFactors = FALSE)
  alignment$genome <- plyr::mapvalues(x = alignment$genome, from = plotting_names$genome, to = plotting_names$plotting_name)
  
  # Truncate
  flog.info("Trimming MSA")
  alignment_trunc <- dplyr::filter(alignment, position <= params$msa_end_pos & position >= params$msa_start_pos)
  alignment_summ_trunc <- dplyr::filter(alignment_summ, position <= params$msa_end_pos & position >= params$msa_start_pos)
  
  # Load residue colours
  res_colours <- read.table("../../residue_colours_rasmol.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                            comment.char = "")
  res_colours$residue <- bio3d::aa321(res_colours$residue_long)
  res_colours <- dplyr::arrange(res_colours, residue)
  
  # Plot the MSA itself
  flog.info("Plotting MSA")
  msa_plot <- ggplot(alignment_trunc, aes(x = position, y = genome)) +
    geom_tile(fill = "white") +
    geom_tile(aes(fill = residue), alpha = 0.7) +
    geom_text(aes(label = residue), size = 3) +
    scale_fill_manual(values = res_colours$colour) +
    guides(fill = FALSE) +
    theme_bw() +
    theme(axis.title = element_blank(),
          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(colour = "#d2d2d2", size = 0.5),
          axis.text.y = element_text(colour = "black"))
    
  # Plot of the summary stats
  flog.info("Plotting MSA stats")
  msa_plot_summ <- ggplot(alignment_summ_trunc, aes(x = position)) +
    geom_bar(aes(weight = similarity), fill = "#333333") +
    geom_text(aes(label = label, colour = heme_binding, y = 1.01), size = 7) +
    guides(colour = FALSE) +
    scale_colour_manual(values = c("black", "blue")) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 9),
          axis.title.y = element_text(size = 9)) +
    ylab("Similarity\n(BLOSUM62)") +
    xlab("Alignment position (bp)") +
    ylim(0,1.3)
  
  # Put everything together and save
  flog.info("Combining into final PDF")
  pdf(file = params$output_filepath, width = 12, height = 5.5)
  egg::ggarrange(tree_plot, msa_plot, ggplot(), msa_plot_summ, nrow = 2, 
                                  heights = c(9,1), widths = c(1,2.5), newpage = TRUE)
  dev.off()
  
  flog.info("Done.")
  
}

if (interactive() == TRUE) {
  ####################################################
  # User variables: #################################
  # See proper descriptions of variables under parse_command_line_input
  # Generate empty list
  params <- list()
  setwd("/home/jmtsuji/Research_General/PhD/04b_Metagenome_resequencing_F2015/10_ATLAS_re_analysis/13_cyc2_comparison/04_phylogeny/extra_genomes/cyc2/")
  params$input_phylogenetic_tree_filepath <- "cyc2_phylogeny_vs5.treefile"
  params$input_msa_filename <- "cyc2_whole_sequence_set_ren1_noLepto_aligned.faa"
  params$output_filepath <- "cyc2_phylogeny_vs5b.pdf"
  params$plotting_names_table_filepath <- "../plotting_names.tsv"
  params$bootstrap_cutoff <- 50
  params$msa_start_pos <- 40
  params$msa_end_pos <- 80
  ####################################################
  
  main(params)
}
