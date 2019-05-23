# ABOUT Figure 02 - Fe/S gene pathways among Chlorobia
Copyright Jackson M. Tsuji, Neufeld Research Group, 2019
Part of the larger *IISD-ELA Chlorobia cyc2 project*.

## 1. Input files
The `input_files` directory contains the data used to create the plot:
- `gene_naming_info.tsv` - gives short abbreviations for all genes to be plotted and the plotting order. Also contains my personal notes.
- `Chlorobia_naming_info.tsv` - links the tree tip names and backblast names to the final names to be plotted
- `chlorobia_riboprotein_phylogeny.treefile` - the phylogenetic tree generated in `Data_analysis_pipeline/06_comparative_genomics/03_chlorobia_phylogeny`
- `backblast_results/*.csv` - the CSV files output from BackBLAST in `Data_analysis_pipeline/06_comparative_genomics/05_pathway_analysis`

## 2. Plotted tree and heatmap using Figure_02_plotter.R
Ran `plot/Figure_02_plotter.R` in interactive mode (e.g., in RStudio) to produce `plot/Figure_02_raw.pdf`. Note that you'll need to install all libraries loaded at the top of the script. After running, I then cleaned up the raw figure in Inkscape to make `plot/Figure_02_cleaned.pdf`, the final figure. Note that the script also outputs `plot/Figure_02_plotting_data.tsv` as a summary of the data to be plotted for the heatmap.

See R package versions in `R_session_info.log`. Log was generated after running the above script by:
```R
sink("R_session_info.log")
sessionInfo()
sink()
```

Note: During iterative testing of the e-value cutoffs, I realized that a few genes had very poor hit profiles on reference organisms. They might be unreliable genes for reciprocal BLAST comparison. I commented them out of the final gene table so that they were excluded from the final plot:
- _qmoA_
- _cysA_
- _cysG_
- _dsrT_
- _soxJ_ and _soxK_ -- N.B., _soxJ_ is a cytochrome
- _dsrJ_
- Also the four proteins with unknown function but possible roles in H2 and sulfite reduction; these are not as relevant to the paper.

