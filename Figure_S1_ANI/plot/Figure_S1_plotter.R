# Nitrospira ANI plotter
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
library(treeio)
library(ggtree, quietly = TRUE, warn.conflicts = FALSE)
library(ape, warn.conflicts = FALSE)
library(maps, warn.conflicts = FALSE)
library(phytools, warn.conflicts = FALSE)
library(reshape2, warn.conflicts = FALSE)
library(gridExtra, warn.conflicts = FALSE)
library(egg, warn.conflicts = FALSE)

################################
# User variables
fastani_data_filepath <- here::here("input_data", "Chlorobia_FastANI_results.txt")
phylogenetic_tree_filepath <- here::here("input_data", "chlorobia_riboprotein_phylogeny.txt")
Chlorobia_naming_table_filepath <- here::here("input_data", "Chlorobia_naming_info.txt")
output_pdf_filepath <- here::here("plot", "Chlorobia_ANI_raw.pdf")
bootstrap_cutoff <- 50
tree_root <- "Ignavibacterium_album_JCM_16511_outgroup"
################################

##### FUNCTIONS #####
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

# Function to read a tibble with some simple presets
read_tibble <- function(table_filepath, header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#") {
  data_table <- read.table(table_filepath, header = header, sep = sep, stringsAsFactors = stringsAsFactors,
                           comment.char = comment.char) %>%
    tibble::as_tibble()
  
  return(data_table)
}

# Function to write a table with some simple presets
write_table <- function(x = table_data, file = table_filepath, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE) {
  write.table(x = x, file = file, sep = sep, quote = quote, row.names = row.names, col.names = col.names)
}
################

# Load ANI table and add headers to ANI data based on tool manual
flog.info(glue::glue("Loading ANI table '", fastani_data_filepath, "'"))
ANI_data <- read_tibble(fastani_data_filepath, header = FALSE)
colnames(ANI_data) <- c("query_genome", "reference_genome", "ANI_value", "n_fragment_mappings", "n_query")

# Load Chlorobia naming table
flog.info(glue::glue("Loading Chlorobia naming table '", Chlorobia_naming_table_filepath, "'"))
naming_data <- read_tibble(Chlorobia_naming_table_filepath)

# Some quick math (didn't end up using this in the final plot, but it is interesting)
ANI_data$percent_fragments_mapped <- ANI_data$n_fragment_mappings / ANI_data$n_query * 100

# Load phylogenetic tree
flog.info(glue::glue("Loading phylogenetic tree '", phylogenetic_tree_filepath, "'"))
phylo_tree <- ggtree::read.tree(phylogenetic_tree_filepath)

# Re-root the tree
#tip_label_index <- match(x = tree_root, table = phylo_tree$tip.label)
#phylo_tree <- ggtree::reroot(phylo_tree, node = tip_label_index)
phylo_tree <- treeio::root(phylo_tree, outgroup = tree_root)

# Plot the tree
flog.info("Generaring the tree plot")
tree_plot <- ggtree(phylo_tree, size = 0.8, colour = "black", ladderize = TRUE,
                    branch.length = 0.1) +
  geom_treescale(x = 0.2, y = 11, linesize = 0.8, fontsize = 3, offset = 0.2, width = 0.08) +
  geom_tiplab(align = TRUE, linetype = "dotted", linesize = 0.35,
              size = 0, offset = 0.1) +
  geom_text2(aes(subset = (grepl(pattern = "^[0-9]+$", x = label) & !(isTip) & as.numeric(label) >= bootstrap_cutoff), 
                 label = as.numeric(label)),
             nudge_x = -0.03, nudge_y = 0.6, size = 2) +
  scale_y_discrete(expand = c(0,0.6)) # to manually make it correspond to the heatmap in y-coordinates

# Change the FastANI names to match the tree names
flog.info("Matching FastANI table to tree tips")
# TODO - this assumes all genomes are in the same directory and have the same extension. Prone to failure.
ANI_data$query_genome <- plyr::mapvalues(ANI_data$query_genome, from = naming_data$fastani_name, to = naming_data$tree_name)
ANI_data$reference_genome <- plyr::mapvalues(ANI_data$reference_genome, from = naming_data$fastani_name, to = naming_data$tree_name)

# Get the order of the tips in the tree plot
tip_order <- dplyr::filter(tree_plot$data, isTip == TRUE)
tip_order <- tip_order[order(tip_order$y, decreasing = TRUE),]$label

# Remove anything in the ANI but not in the tree
if (length(dplyr::setdiff(x = ANI_data$query_genome, y = tip_order)) > 0) {
  ani_entries_not_in_tree <- dplyr::setdiff(x = ANI_data$query_genome, y = tip_order)
  flog.warn(glue::glue("Removing ", length(ani_entries_not_in_tree), " entries from final plot (missing in tree): ", 
                       glue::glue_collapse(ani_entries_not_in_tree, sep = ", ")))
  ANI_data <- dplyr::filter(ANI_data, !(query_genome %in% ani_entries_not_in_tree)) %>%
    dplyr::filter(!(reference_genome %in% ani_entries_not_in_tree))
  naming_data <- dplyr::filter(naming_data, !(subject_name %in% ani_entries_not_in_tree))
}

# Check that the exact same genome names exist between the FastANI data and the tree
if (is.character(all.equal(sort(tip_order), sort(unique(ANI_data$query_genome)))) == TRUE) {
  flog.warn("FastANI labels and tree tips DO NOT MATCH. Results will be unreliable!")
  flog.warn(print(all.equal(sort(tip_order), sort(unique(ANI_data$query_genome)))))
} else if (is.character(all.equal(sort(tip_order), sort(unique(ANI_data$reference_genome)))) == TRUE) {
  flog.warn("FastANI labels and tree tips DO NOT MATCH. Results will be unreliable!")
  flog.warn(print(all.equal(sort(tip_order), sort(unique(ANI_data$reference_genome)))))
}

# Order the genomes in the FastANI table to correspond to the tree tip order
ANI_data$query_genome <- factor(ANI_data$query_genome, levels = rev(tip_order), ordered = TRUE)
ANI_data$reference_genome <- factor(ANI_data$reference_genome, levels = rev(tip_order), ordered = TRUE)

# Change taxon names to those to appear in the final plot
ANI_data$query_genome <- plyr::mapvalues(ANI_data$query_genome, from = naming_data$tree_name, to = naming_data$plotting_name)
ANI_data$reference_genome <- plyr::mapvalues(ANI_data$reference_genome, from = naming_data$tree_name, to = naming_data$plotting_name)


# Add NA values for missing grid values so that grid lines will appear in the final plot
# TODO - A bit hacky
ANI_data_simplified <- reshape2::dcast(ANI_data, query_genome ~ reference_genome, value.var = "ANI_value") %>%
  reshape2::melt(na.rm = FALSE, id.vars = c("query_genome"), variable.name = "reference_genome",
                 value.name = "ANI_value") %>%
  tibble::as_tibble()

# Plot ANI heatmap
flog.info("Plotting the ANI heatmap")
ANI_heatmap <- ggplot2::ggplot(ANI_data_simplified, aes(y = query_genome, x = reference_genome)) +
  geom_tile(aes(fill = ANI_value), colour = "black") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_text(size = 12),
        panel.border = element_rect(colour = "black", size = 1),
        axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.ticks = element_line(size = 0.5), axis.line = element_line(colour = "black", size = 0.5),
        legend.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 10, face = "bold"),
        legend.key = element_rect(colour = "transparent")) +
  guides(fill = guide_legend(title = "Average nucleotide \nidentity (%)")) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(name = "Blues", n = 5), 
                       na.value = "transparent") +
  xlab(NULL) +
  ylab(NULL)

# Put everything together and save
flog.info(glue::glue("Combining tree and heatmap; saving to '", output_pdf_filepath, "'"))
pdf(file = output_pdf_filepath, width = 11, height = 6.5, onefile = FALSE)
egg::ggarrange(tree_plot, ANI_heatmap, nrow = 1,
               widths = c(1,1.5))
dev.off()
