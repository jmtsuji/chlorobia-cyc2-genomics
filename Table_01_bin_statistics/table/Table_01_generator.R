# Table_01_generator.R
# Generates Table 01 for the Chlorobia cyc2 paper
# Copyright Jackson M. Tsuji, 2019

# Load the library and dependencies
library(here)
library(futile.logger)
library(tibble)
library(plyr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(tools)

#################################
## User-defined variables (although more are present in the code below)
statswrapper_filepath <- here::here("input_data", "statswrapper.tsv")
tRNA_counts_filepath <- here::here("input_data", "tRNA_counts.tsv")
checkM_completeness_filepath <- here::here("input_data", "completeness.tsv")
checkM_taxonomy_filepath <- here::here("input_data", "taxonomy.tsv")
chlorobia_naming_guide_filepath <- here::here("input_data", "Chlorobia_naming_guide.tsv")
output_filepath <- here::here("table", "Table_01.tsv")
#################################

# Load and clean up each table
## Statswrapper
flog.info("Loading statswrapper.tsv")
statswrapper_data <- tibble::as_tibble(read.table(statswrapper_filepath, sep = "\t", 
                                                 header = TRUE, stringsAsFactors = FALSE))
statswrapper_data <- dplyr::rename(statswrapper_data, `Bin ID` = filename, Contigs = n_contigs,
                                   `Length (Mb)` = contig_bp, `N50 (kb)` = ctg_L50, L50 = ctg_N50,
                                   `GC content (%)` = gc_avg)
statswrapper_data <- dplyr::select(statswrapper_data, `Bin ID`, Contigs, `Length (Mb)`, `N50 (kb)`,
                                   L50, `GC content (%)`)
statswrapper_data$`Bin ID` <- tools::file_path_sans_ext(statswrapper_data$`Bin ID`)
statswrapper_data$`Length (Mb)` <- round(statswrapper_data$`Length (Mb)` / 10^6, digits = 2)
statswrapper_data$`N50 (kb)` <- round(statswrapper_data$`N50 (kb)` / 10^3, digits = 1)
statswrapper_data$`GC content (%)` <- round(statswrapper_data$`GC content (%)` * 100, digits = 1)

## tRNA counts
flog.info("Loading tRNA_counts.tsv")
tRNA_counts <- tibble::as_tibble(read.table(tRNA_counts_filepath, sep = "\t", 
                                            header = TRUE, stringsAsFactors = FALSE))
tRNA_counts <- dplyr::rename(tRNA_counts, `Bin ID` = Bin_ID)

## CheckM stats
flog.info("Loading completeness.tsv")
checkM_completeness <- tibble::as_tibble(read.table(checkM_completeness_filepath, sep = "\t", 
                                             header = TRUE, stringsAsFactors = FALSE,
                                             comment.char = ""))
checkM_completeness <- dplyr::rename(checkM_completeness, `Bin ID` = "Bin.Id", `Completeness (%)` = Completeness,
                              `Contamination (%)` = Contamination) %>%
                       dplyr::select(`Bin ID`, `Completeness (%)`, `Contamination (%)`)

flog.info("Loading taxonomy.tsv")
checkM_taxonomy <- tibble::as_tibble(read.table(checkM_taxonomy_filepath, sep = "\t", 
                                                header = TRUE, stringsAsFactors = FALSE,
                                                comment.char = ""))
checkM_taxonomy <- dplyr::rename(checkM_taxonomy, `Bin ID` = "Bin.Id", Genes = Gene.count) %>%
                   dplyr::select(`Bin ID`, Genes)

# Join tables
flog.info("Joining tables")
genome_stats <- dplyr::left_join(statswrapper_data, checkM_taxonomy, by = "Bin ID") %>%
                dplyr::left_join(tRNA_counts, by = "Bin ID") %>%
                dplyr::left_join(checkM_completeness, by = "Bin ID")

# Rename genome bins
flog.info("Updating genome names and ordering")
chlorobia_naming_guide <- tibble::as_tibble(read.table(chlorobia_naming_guide_filepath, sep = "\t",
                                                       header = TRUE, stringsAsFactors = FALSE))
genome_stats$`Bin ID` <- plyr::mapvalues(dplyr::pull(genome_stats, `Bin ID`), 
                                         from = dplyr::pull(chlorobia_naming_guide, raw_names),
                                         to = dplyr::pull(chlorobia_naming_guide, final_names))

# Order genome bins
genome_stats$`Bin ID` <- factor(dplyr::pull(genome_stats, `Bin ID`), 
                                levels = dplyr::pull(chlorobia_naming_guide, final_names),
                                ordered = TRUE)
genome_stats <- dplyr::arrange(genome_stats, `Bin ID`)

# Write the table
flog.info("Writing output file")
write.table(genome_stats, output_filepath, sep = "\t", quote = FALSE, 
            col.names = TRUE, row.names = FALSE)
