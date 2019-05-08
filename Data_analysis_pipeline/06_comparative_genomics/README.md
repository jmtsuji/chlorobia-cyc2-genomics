# ABOUT comparative genomics of *Chlorobia* genome bins
Copyright Jackson M. Tsuji, Neufeld Research Group, 2019
Part of the larger *IISD-ELA Chlorobia cyc2 project*.

All code here is to be run in a Bash terminal of a unix system (e.g., Ubuntu 16 or 18) unless indicated otherwise.

Most of the analyses done here generate output files used to generate publication figures. The figure that each output data file is associated with is indicated.

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


## Gene neighbourhood of *Chlorobia* cyc2
Used for Figure 1, panel B


## Phylogeny of *cyc2* genes
Used for Figure 1, panel C


## Phylogenetic placement of *Chlorobia* genome bins
Used for Figure 2
### Download *Chlorobia* genomes


### Identify shared ribosomal proteins


### Generate ribosomal protein phylogeny


## *Chlorobia* gene pathway analysis
Used for Figure 2 heatmap

## Comparative phylogeny of *cyc2* and ribosomal proteins
Used for Supplementary Figure S4



### Done!
This is the end of the main heavy-lifting data processing work for this paper. Figures were generated based off this dataset using code found in each figure folder.

