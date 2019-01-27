# ABOUT assembly and binning of lake metagenome samples
Copyright Jackson M. Tsuji, Neufeld Research Group, 2019
Part of the larger *IISD-ELA Chlorobia cyc2 project*.

All code here is to be run in a Bash terminal of a unix system (e.g., Ubuntu 16 or 18) unless indicated otherwise.

## 1. Data access
Sequencing was performed on an Illumina HiSeq 2500 with 2x200 bp reads (see Methods). Raw reads for the project can be downloaded from NCBI using the `` tool:

```

```

## 2. Assembly via the ATLAS pipeline
To install ATLAS 1.0.22:
```

```

Supplementary file ___ is the config file used for the ATLAS run. You will need to edit the filepaths to the samples to match their locations on your system.

To start the ATLAS run:
```

```
__WARNING__: You'll need a server with at least ~100 GB of RAM to do this. The run might take about a week to finish.

Alternatively, the contigs from the ATLAS run are available on NCBI by running the following:
```

```

## 3. Co-assembly and differential abundance binning of selected samples
A wrapper around the ATLAS pipeline, [co-assembly.sh](___) (version ___) was used to perform co-assembly and is runnable via Docker container. The workflow is described in the Methods section of the paper. To start the Docker container:
```

```

To start the run:
```

```
__WARNING__: You'll need a server with at least ~100 GB of RAM to do this. The run might take about a week to finish.

Alternatively, the contigs from the ATLAS run are available on NCBI by running the following:
```

```

## 4. Bin cleanup
Genome bins were imported into Anvi'o for manual cleanup using [atlas-to-anvi.sh](___), version ____.

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

## 5. Bin dereplication
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

The output table here is included in `04_extra_files` as `____` for interest. It was cleaned up to produce Supplementary File ___, with the addition of the enrichment culture bins.

## Next steps
All dereplicated high-quality _Chlorobia_ bins were manually cleaned up using Anvi'o.
See the tables and figures folder for cyc2-specific analyses steps taken.


