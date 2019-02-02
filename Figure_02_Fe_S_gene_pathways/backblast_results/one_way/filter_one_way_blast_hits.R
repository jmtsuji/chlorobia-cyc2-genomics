# Chlorobia cyc2 manuscript, Figure 02 - support script
# Copyright Jackson M. Tsuji, 2019
# Neufeld Research Group
# Filters one-way BLAST hits for reference organisms.

# Load libraries
library(tibble)
library(dplyr)
library(futile.logger)
library(glue)

######################
## User-defined variables
identity <- 20
query_subject_filename_separator <- "__to__"
input_table_filepaths <- list.files(path = here::here("backblast_results", "one_way"), pattern = "*.csv", full.names = TRUE)
######################

# Function to filter the tables by pident and by best hit
filter_blast_hits <- function(table_filepath) {
  data_table <- tibble::as_tibble(read.table(table_filepath, sep = ",", header = FALSE, stringsAsFactors = FALSE))
  colnames(data_table) <- c("qseqid", "sseqid", "pident", "evalue", "qcovhsp", "bitscore")
  
  # Get subject and query names
  table_filename_base <- basename(tools::file_path_sans_ext(table_filepath))
  table_filename_split <- strsplit(table_filename_base, split = query_subject_filename_separator)[[1]]
  if (length(table_filename_split) != 2) {
    flog.fatal(glue::glue(basename(table_filepath), ": splitting filename by '", query_subject_filename_separator, 
                          "' results in ", length(table_filename_split), " entries (looking for 2). ",
                          "Cannot parse query/subject names. Exiting..."))
    quit(status = 1)
  }
  query_name <- table_filename_split[1]
  subject_name <- table_filename_split[2]
  
  # No need to add special columns to the table for the subject and query, because it will be saved in the standaard BackBLAST output format that lacks this info
  
  # Filter by percent identity cutoff
  data_table <- dplyr::filter(data_table, pident > identity)
  
  # Filter to best hit
  data_table <- dplyr::arrange(data_table, evalue)
  data_table <- dplyr::distinct(data_table, qseqid, .keep_all = TRUE)
  
  # Warn user if any did not hit to self -- will need to check these manually
  for (i in 1:nrow(data_table)) {
    if (data_table$qseqid[i] != data_table$sseqid[i]) {
      flog.warn(glue::glue("qseqid '", data_table$qseqid[i], "' did not have best hit to self. Check on this entry manually."))
    }
  }
  
  flog.info("Saving filtered CSV file")
  write.table(data_table, file = here::here("backblast_results", "one_way", paste(table_filename_base, ".csv.filtered", sep = "")), 
              sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  return(data_table)
}

# Run on both tables. Auto-saves the output files.
data_tables <- lapply(input_table_filepaths, filter_blast_hits)
