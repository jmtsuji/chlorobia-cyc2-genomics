# Figure_2A_plotter.R
# Generates panel A of Supplementary Figure 2 for the Chlorobia cyc2 paper
# Copyright Jackson M. Tsuji, 2020
# N.B., Contains hard-coded tweaks for c5 the family cytochrome alignment; don't use arbitrarily for any MSA

# Load required packages
library(here)
library(futile.logger)
library(plyr)
library(tibble)
library(dplyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(gridExtra, warn.conflicts = FALSE)
library(egg, warn.conflicts = FALSE)
library(bio3d)

#####################################################
## User-defined variables
#####################################################
KEY_RESIDUES <- c(95,98,99)

params <- list()

params$input_msa_filename <- here::here("input_data", "c5_family_Chlorobia_aligned.faa")
params$output_filepath <- here::here("plot", "Figure_S2A_raw.pdf")
params$plotting_names_table_filepath <- here::here("input_data", "c5_plotting_name_info.tsv")
params$msa_residue_colour_guide_filepath <- here::here("input_data", "residue_colours_rasmol.tsv")
#####################################################

# Function to convert an msa object from bio3d into a data frame that works well with ggplot
create_alignment_data_frame <- function(msa_data) {
  # Convert to data frame
  alignment <- tibble::as_tibble(t(msa_data$ali))
  alignment$position <- 1:nrow(alignment)
  alignment <- tibble::as_tibble(reshape2::melt(alignment, id.vars = "position", variable.name = "genome", value.name = "residue"))
  alignment$genome <- as.character(alignment$genome)
  alignment <- dplyr::filter(alignment, residue != "-") # Remove dashes
  
  return(alignment)
}

# Function to calculate similarity of the alignment and annotate key residues
# Some HARD-CODED numbers are in here, N.B.
generate_alignment_summary_stats <- function(alignment, msa_data) {
  # Make summary
  alignment_summ <- tibble::tibble(position = 1:max(alignment$position))
  alignment_summ$similarity <- bio3d::conserv(msa_data, method = "similarity", sub.matrix = "blosum62")
  alignment_summ$label <- unlist(lapply(alignment_summ$similarity, function(x) {
    if (x == 1) {
      return("*")
    } else {
      return(NA)
    }
  }))
  # Manually add predicted H25 location
  alignment_summ$key_residue <- logical(length = nrow(alignment_summ))
  alignment_summ$key_residue[KEY_RESIDUES] <- "TRUE"
  alignment_summ$label[KEY_RESIDUES] <- "*"
  
  return(alignment_summ)
}

main <- function(params) {
  # Load MSA data and process
  flog.info("Loading MSA")
  msa_data <- bio3d::read.fasta(file = params$input_msa_filename, rm.dup = TRUE, to.dash = TRUE)
  alignment <- create_alignment_data_frame(msa_data)
  alignment_summ <- generate_alignment_summary_stats(alignment, msa_data)
  
  # Also change names to more readable
  plotting_names <- tibble::as_tibble(read.table(params$plotting_names_table_filepath, sep = "\t", 
                                                 header = TRUE, stringsAsFactors = FALSE))
  alignment$genome <- plyr::mapvalues(x = alignment$genome, from = plotting_names$genome, 
                                      to = plotting_names$plotting_name)
  
  # Define tip order
  alignment$genome <- factor(alignment$genome, levels = rev(plotting_names$plotting_name), ordered = TRUE)
  
  # Load residue colours
  res_colours <- tibble::as_tibble(read.table(params$msa_residue_colour_guide_filepath, sep = "\t", 
                                              header = TRUE, stringsAsFactors = FALSE, comment.char = ""))
  # Get the matching 1-letter code for each three-letter abbreviation
  res_colours$residue <- bio3d::aa321(res_colours$residue_long)
  res_colours <- dplyr::arrange(res_colours, residue)
  
  # Plot the MSA itself
  flog.info("Plotting MSA")
  msa_plot <- ggplot(alignment, aes(x = position, y = genome)) +
    geom_tile(fill = "white") +
    geom_tile(aes(fill = residue), alpha = 0.7) +
    #geom_text(aes(label = residue), size = 3, family = "Courier") +
    scale_fill_manual(values = res_colours$colour) +
    guides(fill = FALSE) +
    theme_bw() +
    theme(axis.title = element_blank(),
          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(colour = "#d2d2d2", size = 0.5),
          axis.text.y = element_text(colour = "black"))
  
  # Plot the summary stats
  flog.info("Plotting MSA stats")
  msa_plot_summ <- ggplot(alignment_summ, aes(x = position)) +
    geom_bar(aes(weight = similarity), fill = "#333333") +
    geom_text(aes(label = label, colour = key_residue, y = 1.01), size = 7) +
    guides(colour = FALSE) +
    scale_colour_manual(values = c("black", "red")) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 9),
          axis.title.y = element_text(size = 9)) +
    ylab("Similarity\n(BLOSUM62)") +
    xlab("Alignment position (amino acid)") +
    ylim(0,1.5)
  
  # Put everything together and save
  flog.info("Combining into final PDF")
  pdf(file = params$output_filepath, width = 8, height = 2.8, onefile = FALSE)
  egg::ggarrange(msa_plot, msa_plot_summ, nrow = 2, 
                 heights = c(4,1), newpage = TRUE)
  dev.off()
  
  flog.info("Done.")
  
}

main(params)
