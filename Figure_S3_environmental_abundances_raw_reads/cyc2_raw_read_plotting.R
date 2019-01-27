# Plotting MetAnnotate raw read data for the Chlorobia cyc2 paper (supplementary figure)
# Copyright Jackson M. Tsuji, 2019
# Relies on the metannotate-analysis library, version 0.9.5 (see link below)

# Load the library and dependencies
source("https://github.com/jmtsuji/metannotate-analysis/releases/download/v0.9.5/metannotate_barplots.R")

#################################
## User-defined variables (although more are present in the code below)
setwd("/home/jmtsuji/Research_General/PhD/04b_Metagenome_resequencing_F2015/10_ATLAS_re_analysis/07b_metannotate_raw_reads/v0.9.5/")
metannotate_table_filename <- "all_annotations_w9HU0A510082988.tsv.gz"
hmm_naming_info_filename <- "hmm_info_vs7.tsv"
dataset_naming_info_filename <- "dataset_info_vs1.tsv"
output_filename <- "Chlorobia_cyc2_raw_read_plot_vs1_evalue_1e-10_1perc.pdf"
evalue <- 1e-10
plot_width <- 200
plot_height <- 250
#################################

# Load the MetAnnotate table
metannotate_data <- read_metannotate_tibble(metannotate_table_filename)

# Generate templates to be filled in by the user
# N.B., Only really need to run once!
setup_templates <- create_setup_templates(metannotate_data, write_tables = TRUE)
# Saves as 'hmm_info_template.tsv' and 'dataset_info_template.tsv'.
# I edit those to be the 'hmm_naming_info_filename' and 'dataset_naming_info_filename' above.

# Map the edited template info onto the MetAnnotate table
metannotate_data_mapped <- map_naming_information(metannotate_data, hmm_naming_info_filename, dataset_naming_info_filename)

# Generate the bubble plot (N.B., LOTS of hard-coded settings here)
metannotate_plot <- explore_metannotate_data(metannotate_data_mapped, evalue = evalue, taxon = "Family",
                                             normalizing_HMM = "rpoB", plot_type = "bubble", 
                                             top_x = 0.01, percent_mode = "within_sample",
                                             colouring_template_filename = NA,
                                             space = "free", bubble_size_range = c(1,50), alpha = 0.5,
                                             bubble_labels = TRUE)
print(metannotate_plot)


### I'm choosing not to specify custom plotting colours for this workflow.
# metannotate_plot <- explore_metannotate_data(metannotate_data_mapped, evalue = 1e-10, taxon = "Family",
#                                              normalizing_HMM = "rpoB", plot_type = "bar", 
#                                              top_x = 10, percent_mode = "within_sample",
#                                              colouring_template_filename = "colouring_template.tsv")


# metannotate_plot <- explore_metannotate_data(metannotate_data_mapped, evalue = 1e-10, taxon = "Family",
#                                              normalizing_HMM = "rpoB", plot_type = "bar", 
#                                              top_x = 10, percent_mode = "within_sample",
#                                              colouring_template_filename = "colouring_guide.tsv")
# print(metannotate_plot)

# Save the plot
ggsave(file = output_filename, width = plot_width, 
       height = plot_height, units = "mm")

