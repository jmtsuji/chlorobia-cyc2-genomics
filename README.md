# ABOUT IISD-ELA _Chlorobia_ _cyc2_ genomics project
Copyright Jackson M. Tsuji, Neufeld Research Group, 2020  
[![Github repo DOI](https://zenodo.org/badge/168807016.svg)](https://zenodo.org/badge/latestdoi/168807016)

This Github repository describes the steps/code used to perform bioinformatics analysis for the paper by Tsuji and colleagues, "Genomic potential for photoferrotrophy in a seasonally anoxic Boreal Shield lake" (2019; [biorxiv pre-print available here](https://doi.org/10.1101/653014). Descriptor files in this repo are in the Markdown (`.md`) format and begin with `README`. Other files include scripts and raw data files.

## The information is organized as follows:
### `Data_analysis_pipeline`
Describes the protocol used to produce assembled contigs and genome bins from the raw metagenomic datasets used for this paper.

Includes (some highlights):
- `01_data_acquisition` - contains instructions for how to download the online sequencing data associated with this paper, including raw read data and genome bins
- `04_HMM_development/02_alignment` - contains an alignment of primary sequences of the *Chlorobia*-associated *cyc2* genes identified in this study, plus reference *cyc2* sequences
- `04_HMM_development/03_HMMs` - contains the custom HMMs developed in this study
- `05_bin_analysis/03_metannotate/genome_bins` - contains primary sequences of non-*Chlorobia* *cyc2* genes detected in the lake metagenomes

### Tables and figures
A unique folder exists for most main or supplementary tables/figures in the publication, e.g., `Figure_01_cyc2_genes`. These describe the steps required to analyze the data for the figure and to generate a raw version of the figure for downstream editing.

### `Other`
Some supplementary material did not warrant having its own folder, so it is grouped in the `Other` folder in the same way that it would appear in the Supplementary Information on the journal website.

## Additional information about this repo
### Zenodo
Note that this repo has a corresponding data repository in Zenodo for hosting larger files (e.g., non-curated genome bins). Check out the Zenodo version (`v1.0.1`) corresponding to the current Gitub repo code below:  
[![Zenodo data repository](https://zenodo.org/badge/DOI/10.5281/zenodo.3228469.svg)](http://doi.org/10.5281/zenodo.3228469)

Code in this repo shows how to use the files in the Zenodo repo for analysis.

### Issues
It's likely that there are some typo's in the code. I did not perform all of these analyses at once, but I have tried in making this repo to organize them into one smooth, continuous workflow. This involved renaming some variables and such. Not all of the code has been tested in its current state.

If you want to try to use something here and are having problems (or have other questions/concerns about the analysis), feel free to post a [Github issue](https://github.com/jmtsuji/chlorobia-cyc2-genomics/issues).

### Citation
If you'd like to cite this repo, then please cite the [biorxiv preprint](https://doi.org/10.1101/653014).

## Final remarks
I hope this is a helpful resource for the scientific community, both in the critique of my work and in further development of bioinformatics. Feel free to post an issue if you have any questions or concerns. Enjoy!  


Jackson M. Tsuji  
PhD Candidate  
Neufeld Research Group  
University of Waterloo  
Waterloo, Ontario, CANADA
