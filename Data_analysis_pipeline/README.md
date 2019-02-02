# ABOUT assembly and binning of lake metagenome samples
Copyright Jackson M. Tsuji, Neufeld Research Group, 2019
Part of the larger *IISD-ELA Chlorobia cyc2 project*.

## How this folder works
This folder of the project contains the code that I used to do most of the heavy-lifting of metagenomic data analysis, along with a few support files. If you follow along with the code as given in each sub-module's `README.md` file, you should be able to reproduce my results.

## How do I use this folder?
Clone the whole repo onto your Linux server, change the `zenodo_repo_location` variable at the top of each submodule's `README.md` file, and then run the code as-is. Each step in the pipeline will be executed within the repo folder. Running the code here will download the metagenomes (`01_data_acquisition`), assemble/bin them (`02_assembly_and_binning`), and dereplicate them and import them for curation (`03_bin_curation`).

## Do I want to run your code on my server?
Probably not. The whole process requires at least ~100 GB of RAM and at least a week of compute time. The main reason for me including the code here is for readers to understand the workflow that I used, not necessarily to re-run it themselves. __In fact, I expect much of the code here to become obsolete within the next few weeks/months as ATLAS and other metagenomics pipelines continue to improve.__ So even if you want to run similar analyses to mine, I would recommend trying (ATLAS version 2)[https://github.com/metagenome-atlas/atlas] as soon as it is available rather than reproducing my workflow.

## Requirements to run the code
- All code here is to be run in a Bash terminal of a unix system (e.g., Ubuntu 16 or 18) unless indicated otherwise
- You'll need at least ~300 GB of storage, ~100 GB of RAM, and possibly a week or two of compute time (e.g., on 12 threads)
- The only depedencies are `conda` (recommended is `miniconda3`) and `docker`
- I've tried to write this documentation to be accessible to an 'early intermediate'-level Bash coder with knowledge of metagenomics

## To just get the data
You might want to skip all this and just access the data itself. I've uploaded data from several steps along the pipeline for ease of public use.

### A. The raw reads
Follow the code in the `01_data_acquisition` folder to get the raw metagenomic reads (Illumina).

### B. The assembled metagenomes
Assembled metagenomes have been deposited at ___TODO___

```

```

### C. Curated Chlorobia bins
Curated _Chlorobia_ bins have been deposited at ___TODO___

```

```

### D. All the bins
I did not upload all of the bins onto NCBI because they have not been curated and may contain contaminants. However, they are available in the (Zenodo repository at this link)[]. Download them by:

```

```

## Conclusions
I hope this is helpful to you! Please see the other folders is this repo for code used to analyze the metagenomes/bins and figure generation.

