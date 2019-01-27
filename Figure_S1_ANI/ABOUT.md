# ABOUT Supplementary Figure S1 - average nucleotide identity among Chlorobia
Copyright Jackson M. Tsuji, Neufeld Research Group, 2019
Part of the larger *IISD-ELA Chlorobia cyc2 project*.

All code here is to be run in a Bash terminal of a unix system (e.g., Ubuntu 16 or 18) unless indicated otherwise.

# Choose your working directory -- you need to modify this variable based on where you download the Zenodo repository to:
For example:
```
work_dir="/Analysis/jmtsuji/Chlorobia_cyc2_code/Figure_S1_ANI"
cd ${work_dir}
```

## Part A: FastANI analysis
All of Part A is run in the Linux command line (Ubuntu 18.04) via Bash. The only dependency is [miniconda](https://docs.conda.io/en/latest/miniconda.html).

### 1. Created the conda environment and made directory structure
```
conda create -n fastani -c bioconda fastani=1.1
cd ${work_dir}
mkdir -p genome_data plotting_data
```

### 2. Collected genomes
N.B., the genome_data folder has been left empty in the repo to save space. You need to run this code to re-download.
```
# Same code as for Figure 2

# TODO


### Reference genomes


### Genome bins
```

### 3. Ran FastANI
```
# Made reference and query lists for FastANI; both are the same because I am doing an all-by-all comparison
realpath genome_data/*fa > plotting_data/query.list
cp plotting_data/query.list plotting_data/reference.list

# Ran FastANI
fastANI --rl plotting_data/reference.list --ql plotting_data/query.list -o plotting_data/Chlorobia_FastANI_results.txt 2>&1 | tee plotting_data/Chlorobia_FastANI.log
# Completes in less than a minute, even running on a single thread.

# Removed temp files
rm plotting_data/reference.list plotting_data/query.list
```

## Part B: Visualized in R
The following three input data files are required:
- FastANI output from above
- Phylogenetic tree of _Chlorobia_ based on concatenated ribosomal proteins -- same as made for Figure 2, i.e., `Chlorobia_riboprotein_tree.treefile`
- Made a table of input genome filepaths and genome names (matching phylogenetic tree tip labels). Saved as `Chlorobia_naming_info.tsv`. First column (`fastani_name`) is the full filepath to the genome file. Second column (`tree_name`) is the genome name in the phylogenetic tree. Third column (`plotting_name`) is the final name to appear on the plot.
All three data files are saved in the `plotting_data` sub-directory.

Running the R script `Figure_S1_plotter.R` in interactive mode (e.g., in RStudio) will generate the base figure `Figure_S1_raw.pdf`. Note that you'll also need to install all libraries loaded at the top of the script (e.g., via `install.packages()` and [`BiocManager::install()`](https://bioconductor.org/install/) in the R console).
I cleaned up the raw plot in Inkscape to make `Figure_S1_cleaned.pdf`, the final version that appears in the paper.

