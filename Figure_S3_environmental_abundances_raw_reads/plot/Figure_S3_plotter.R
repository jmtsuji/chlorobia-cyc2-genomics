# Plotting MetAnnotate raw read data for the Chlorobia cyc2 paper (supplementary figure S3)
# Copyright Jackson M. Tsuji, 2019
# Relies on the metannotate-analysis library, version 0.9.5 (see link below)

# Load the library and dependencies
library(here)
source("https://github.com/jmtsuji/metannotate-analysis/releases/download/v0.9.5/metannotate_barplots.R")

#################################
## User-defined variables (although more are present in the code below)
metannotate_table_filename <- here::here("input_data", "metannotate_annotations_unassembled_reads.tsv.gz")
hmm_naming_info_filename <- here::here("input_data", "hmm_info.tsv")
dataset_naming_info_filename <- here::here("input_data", "dataset_info.tsv")
output_basename <- here::here("plot", "Figure_S3")
abundance_threshold <- 0.01 # Expressed as a proportion; any taxon under this will not be plotted
plotting_taxon <- "Family"
evalue <- 1e-10
plot_width <- 200
plot_height <- 250
print_templates <- FALSE # I already printed guide templates and modified them to make 'hmm_naming_info_filename' and 'dataset_naming_info_filename'
#################################

# Load the MetAnnotate table
metannotate_data <- read_metannotate_tibble(metannotate_table_filename)

# Generate templates to be filled in by the user
if (print_templates == TRUE) {
  flog.info("Printing setup templates")
  setup_templates <- create_setup_templates(metannotate_data, write_tables = TRUE)
  flog.info("Setup templates printed to 'hmm_naming_info_filename' and 'dataset_naming_info_filename'. Please edit these and then use to run this script. Exiting...")
  quit(save = "no", status = 0)
}

# Map the edited template info onto the MetAnnotate table
metannotate_data_mapped <- map_naming_information(metannotate_data, hmm_naming_info_filename, dataset_naming_info_filename)

# Generate the bubble plot (N.B., some hard-coded settings here)
metannotate_plot <- explore_metannotate_data(metannotate_data_mapped, evalue = evalue, taxon = plotting_taxon,
                                             normalizing_HMM = "rpoB", plot_type = "bubble", 
                                             top_x = abundance_threshold, percent_mode = "within_sample",
                                             colouring_template_filename = NA,
                                             space = "free", bubble_size_range = c(1,50), alpha = 0.5,
                                             bubble_labels = TRUE)
# Save the plot
ggsave(file = paste(output_basename, ".pdf", sep = ""), width = plot_width, 
       height = plot_height, units = "mm")

# Save a copy of the plotting data for reference
write.table(x = dplyr::select(metannotate_plot$data, -label), 
            file = paste(output_basename, "_plotting_data.tsv", sep = ""),
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

####### Do this again, for reference, with a more relaxed percent abundance threshold
abundance_threshold_relaxed <- 0.005
metannotate_plot_relaxed <- explore_metannotate_data(metannotate_data_mapped, evalue = evalue, taxon = plotting_taxon,
                                             normalizing_HMM = "rpoB", plot_type = "bubble", 
                                             top_x = abundance_threshold_relaxed, percent_mode = "within_sample",
                                             colouring_template_filename = NA,
                                             space = "free", bubble_size_range = c(1,50), alpha = 0.5,
                                             bubble_labels = TRUE)
ggsave(file = paste(output_basename, "_0.5percent_for_reference.pdf", sep = ""), width = plot_width, 
       height = plot_height, units = "mm")
write.table(x = dplyr::select(metannotate_plot_relaxed$data, -label), 
            file = paste(output_basename, "_0.5percent_for_reference_plotting_data.tsv", sep = ""),
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
