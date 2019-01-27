# Chlorobia cyc2 manuscript, Figure S1 plotter
# Copyright Jackson M. Tsuji, 2019
# Neufeld Research Group

# Load libraries
library(here)
library(futile.logger)
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

################################
# HARD-CODED plotting variables
# setwd("/home/jmtsuji/Research_General/PhD/04b_Metagenome_resequencing_F2015/10_ATLAS_re_analysis/15_final_zenodo_code/Chlorobia_cyc2_code/Figure_S1_ANI/")
fastani_data_filepath <- here::here("plotting_data", "Chlorobia_FastANI_results.txt")
phylogenetic_tree_filepath <- here::here("plotting_data", "Chlorobia_riboprotein_tree.treefile")
Chlorobia_naming_table_filepath <- here::here("plotting_data", "Chlorobia_naming_info.tsv")
output_pdf_filepath <- here::here("Figure_S1_raw.pdf")
bootstrap_cutoff <- 50
tree_root <- "Ignavibacterium_album_JCM_16511_outgroup"
################################

# Load ANI table and add headers to ANI data based on tool manual
flog.info(glue::glue("Loading ANI table '", fastani_data_filepath, "'"))
ANI_data <- tibble::as_tibble(read.table(fastani_data_filepath, sep = "\t", header = FALSE, stringsAsFactors = FALSE))
colnames(ANI_data) <- c("query_genome", "reference_genome", "ANI_value", "n_fragment_mappings", "n_query")

# Load Chlorobia naming table
flog.info(glue::glue("Loading Chlorobia naming table '", Chlorobia_naming_table_filepath, "'"))
Chlorobia_names <- tibble::as_tibble(read.table(Chlorobia_naming_table_filepath, sep = "\t", header = TRUE, stringsAsFactors = FALSE))

# Some quick math (didn't end up using this in the final plot, but it is interesting)
ANI_data$percent_fragments_mapped <- ANI_data$n_fragment_mappings / ANI_data$n_query * 100

# Load phylogenetic tree
flog.info(glue::glue("Loading phylogenetic tree '", phylogenetic_tree_filepath, "'"))
phylo_tree <- ggtree::read.tree(phylogenetic_tree_filepath)

# Re-root the tree
tip_label_index <- match(x = tree_root, table = phylo_tree$tip.label)
phylo_tree <- ggtree::reroot(phylo_tree, node = tip_label_index)

# Plot the tree
flog.info("Generaring the tree plot")
tree_plot <- ggtree(phylo_tree, size = 1.5, colour = "black", ladderize = TRUE,
                    branch.length = 0.1) +
  geom_treescale(x = 0.4, y = 10, linesize = 1, fontsize = 3, offset = 0.2, width = 0.2) +
  geom_tiplab(align = TRUE, linetype = "dotted",
              size = 0, offset = 0.1) +
  geom_text2(aes(subset = (grepl(pattern = "^[0-9]+$", x = label) & !(isTip) & as.numeric(label) > bootstrap_cutoff), 
                 label = as.numeric(label)),
             nudge_x = -0.05, nudge_y = 0.4, size = 3) +
  scale_y_discrete(expand = c(0,0.6)) # to manually make it correspond to the heatmap in y-coordinates

# Change the FastANI names to match the tree names
flog.info("Matching FastANI table to tree tips")
ANI_data$query_genome <- plyr::mapvalues(ANI_data$query_genome, from = Chlorobia_names$fastani_name, to = Chlorobia_names$tree_name)
ANI_data$reference_genome <- plyr::mapvalues(ANI_data$reference_genome, from = Chlorobia_names$fastani_name, to = Chlorobia_names$tree_name)

# Get the order of the tips in the tree plot
tip_order <- dplyr::filter(tree_plot$data, isTip == TRUE)
tip_order <- tip_order[order(tip_order$y, decreasing = TRUE),]$label

# Check that the exact same genome names exist between the FastANI data and the tree
if (all.equal(sort(tip_order), sort(unique(ANI_data$query_genome))) == FALSE) {
  flog.warn("FastANI labels and tree tips DO NOT MATCH. Results will be unreliable!")
} else if (all.equal(sort(tip_order), sort(unique(ANI_data$reference_genome))) == FALSE) {
  flog.warn("FastANI labels and tree tips DO NOT MATCH. Results will be unreliable!")
}

# Order the genomes in the FastANI table to correspond to the tree tip order
ANI_data$query_genome <- factor(ANI_data$query_genome, levels = rev(tip_order), ordered = TRUE)
ANI_data$reference_genome <- factor(ANI_data$reference_genome, levels = rev(tip_order), ordered = TRUE)

# Change taxon names to those to appear in the final plot
ANI_data$query_genome <- plyr::mapvalues(ANI_data$query_genome, from = Chlorobia_names$tree_name, to = Chlorobia_names$plotting_name)
ANI_data$reference_genome <- plyr::mapvalues(ANI_data$reference_genome, from = Chlorobia_names$tree_name, to = Chlorobia_names$plotting_name)

# Plot ANI heatmap
flog.info("Plotting the ANI heatmap")
ani_heatmap <- ggplot(ANI_data, aes(x = reference_genome, y = query_genome)) +
  geom_tile(aes(fill = ANI_value)) +
  geom_text(aes(label = round(ANI_value, digits = 1)), size= 2.5) +
  scale_fill_gradient(na.value = "black", low = "#003300", high = "#99ff33") +
  theme_bw() +
  theme(text = element_text(colour = "black", size = 12), line = element_line(colour = "black"),
        rect = element_rect(colour = "black"), panel.grid = element_blank(),
        aspect.ratio = 1,
        axis.text = element_text(size = 10, face = "italic", colour = "black"), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.ticks = element_line(size = 0.5), axis.line = element_line(size = 0.5),
        legend.text = element_text(size = 10), legend.title = element_text(size = 8),
        legend.key = element_rect(colour = "transparent"), legend.key.size = unit(5, "mm")) +
  xlab("") +
  ylab("") +
  guides(fill = guide_legend(title = "Average nucleotide\nidentity (%)"))

# Put everything together and save
flog.info(glue::glue("Combining tree and heatmap; saving to '", output_pdf_filepath, "'"))
pdf(file = output_pdf_filepath, width = 14, height = 7, onefile = FALSE)
egg::ggarrange(tree_plot, ani_heatmap, nrow = 1,
               widths = c(2,2.5))
dev.off()

