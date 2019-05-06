# ABOUT bin dereplication and manual curation
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

```
Use this environment via `` (as shown below).



## Bin dereplication
Unique genome bins were determined using dRep, version 2.0.5.

To install:
```

```

To dereplicate the bins:
```

```

Can run some custom R code afterward to nicely summarize the results, in an R console:
```

```

The output table here is included in `Table_S1.csv`, a copy of which is included in this folder. It was cleaned up to produce Supplementary File ___, with the addition of the enrichment culture bins.


## Manual cleanup of selected bins
All genome bins were imported into Anvi'o for subsequent manual cleanup of specific bins of interest using [atlas-to-anvi.sh](___), version ____.

To install:
```

```

To import the bins:
```

```

To analyze a bin of interest (e.g., `___`):
```

```
Bins of interest were seleted after the bin dereplication step below.


### Done! This is the end of the main heavy-lifting data processing work for this paper. Figures were generated based off this dataset using code found in each figure folder.

