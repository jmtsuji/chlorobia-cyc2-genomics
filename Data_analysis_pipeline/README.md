# ABOUT assembly and binning of lake metagenome samples
Copyright Jackson M. Tsuji, Neufeld Research Group, 2019  
Part of the larger *IISD-ELA Chlorobia cyc2 project*.

## How this folder works
This folder contains the code that I used to do most of the heavy-lifting of metagenomic data analysis, along with a few support files. If you follow along with the code as given in each sub-module's `README.md` file, you should be able to reproduce my results.

The folder also contains instructions for how to download the raw metagenomic reads or genome bins from the project (see `01_data_acquisition`).

## How do I use this folder?
To just download the raw reads: run the code in the `README.md` file in `01_data_acquisition`.

To run the assembly/binning pipeline: clone the whole repo onto your Linux server, change the `github_repo_location` variable at the top of each submodule's `README.md` file, and then run the code as-is. Each step in the pipeline will be executed within the repo folder. Running the code end-to-end will download the metagenomes (`01_data_acquisition`), assemble/bin them (`02_assembly_and_binning`), and dereplicate them and import them for curation (`03_bin_curation`).

## Intended usage
I made my assembly/binning code available so that users/readers and understand and check my code, not so that users can copy this code for their own analyses or re-run my analysis end-to-end. Some reasoning:

- I expect much of this workflow to become obsolete in the coming months. Much of what the code I wrote accomplishes is already built into the latest version of the ATLAS pipeline (2.0.6; I was working with 1.0.22 for this publication). I would recommend trying (ATLAS2)[https://github.com/metagenome-atlas/atlas] rather than modifying my workflow for your data.
- Re-running my analyses end-to-end would be interesting, but keep in mind that it will take at least ~100 GB of RAM possibly a week or more of compute time!

## Requirements to run the code end to end
- All code here is to be run in a Bash terminal of a unix system (e.g., Ubuntu 16 or 18) unless indicated otherwise
- You'll need at least ~300 GB of storage, ~100 GB of RAM, and possibly a week or two of compute time (e.g., on 12 threads)
- The only depedencies are `conda` (recommended is `miniconda3`) and `docker`
- I've tried to write this documentation to be accessible to an 'early intermediate'-level Bash coder with knowledge of metagenomics

## Conclusions
I hope this is helpful to you! Please see the other folders is this repo for code used to analyze the metagenomes/bins and figure generation.

