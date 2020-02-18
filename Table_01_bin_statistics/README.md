# ABOUT Table 1 - genome bin statistics
Copyright Jackson M. Tsuji, Neufeld Research Group, 2019  
Part of the larger *IISD-ELA Chlorobia cyc2 project*.

Raw data tables were compiled to make the final publication table.

## 1. Data collection
The `input_files` directory contains the data used to create the plot:
- `statswrapper.tsv` - basic statistics about each genome generated in `Data_analysis_pipeline/03_bin_curation/04_bin_stats`
- `completeness.tsv` and `taxonomy.tsv` - checkM completeness and taxonomy statistics generated in `Data_analysis_pipeline/03_bin_curation/04_bin_stats`
- `tRNA_counts.tsv` - prokka-based tRNA gene counts for each genome generated in `Data_analysis_pipeline/03_bin_curation/04_bin_stats`
- `Chlorobia_naming_guide.tsv` - guide file specifying the final name of each genome and the desired order in the final table

## 2. Creating the table
Ran `plot/Table_01_generator.R` in interactive mode (e.g., in RStudio) to produce `table/Table_01.tsv`. Note that you'll need to install all libraries loaded at the top of the script. Formatting was cleaned up to make the final Excel version of the table, `table/Table_01.xlsx`. Also note that the name of the Lake 304 enrichment was manually changed to "Ca. Chl. canadense" due to this enrichment culture later being assigned a name.

See R package versions in `R_session_info.log`. Log was generated after running the above script by:
```R
sink("R_session_info.log")
sessionInfo()
sink()
```

