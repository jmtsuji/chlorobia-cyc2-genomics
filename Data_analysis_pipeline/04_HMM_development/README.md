# ABOUT *cyc2* profile Hidden Markov Model development
Copyright Jackson M. Tsuji, Neufeld Research Group, 2019
Part of the larger *IISD-ELA Chlorobia cyc2 project*.

All code here is to be run in a Bash terminal of a unix system (e.g., Ubuntu 16 or 18) unless indicated otherwise.

## Define where you downloaded the Github repo:
```
github_repo_location="/Analysis/jmtsuji/Chlorobia_cyc2_code"
```

## Software prerequisites
- miniconda (miniconda3 preferred)


## Conda environment with all needed dependencies:
```
conda create -y -n dRep -c bioconda -c r drep=2.0.5 r-dplyr
```
Use this environment via `conda activate dRep` (as shown below).


## Download genomes and identify *cyc2* genes


## Align the genes

## Create the HMMs
Ran the classifier (needs ~100 GB RAM!)
```
# Set user variables
input_dir="${github_repo_location}/Data_analysis_pipeline/03_bin_curation/04_all_bin_summary"
output_dir="${github_repo_location}/Data_analysis_pipeline/03_bin_curation/05_gtdbtk_taxonomy"
threads=12

gtdbtk classify_wf --genome_dir ${input_dir} --out_dir ${output_dir} -x fa --min_perc_aa 0 --prefix ELA111314_dRep_gtdbtk --cpus ${threads}
```
The output files `ELA111314_dRep_gtdbtk.ar122.summary.tsv` and `ELA111314_dRep_gtdbtk.bac120.summary.tsv` will be used downstream in Figure 3.

Note: in the actual analysis, GTDBTk was run separately for the lake vs. enrichment culture metagenomes, and the output files were then concatenated. However, running them together is simpler and I think should not change the results much at all.


