# ABOUT Figure 2 - sulfur vs. iron gene pathways in _Chlorobia_
Copyright Jackson M. Tsuji, Neufeld Research Group, 2019
Part of the larger *IISD-ELA Chlorobia cyc2 project*.

The phylogeny and heatmap were produced separately and were manually combined to produce the final figure.

## 1. Data collection
Predicted open reading frames (ORFs) of reference strains in the _Chlorobiaceae_ family, along with an outgroup (_Ignavibacterium album_) were downloaded from NCBI:
```

```

Predicted open reading frames (ORFs) of dereplicated high-quality _Chlorobia_ genome bins from this study (see Assembly and Binning ABOUT files) were also used.

The final set of ORF (`.faa`) files is available in `04_extra_files` in this repository.

## 2. Constructing the concatenated ribosomal protein phylogeny
### Part A - Identifying marker genes from the rp1 ribosomal protein set
This was a semi-manual process.


The phylogeny was visualized using Dendroscope and saved as a PDF. The PDF was then opened using Inkscape to increase line thicknesses, simplify labels, and so on.

## 3. Building the reciprocal BLAST heatmap
### Part A - choosing the gene set and e-value cutoff


### Part B - performing BLAST


### Part C - building the heatmap
Input files required (included in `04_extra_files` in this repository):
- BLAST hit table (___; from Part B above)
- Gene naming/ordering table (___)
- Genome naming/ordering table (___; ordered in the same order as the phylogenetic tree labels)

The following R code was used to build the final heatmap. You'll need to adjust the working directory to match whatever folder you work in. You'll also need to install all libraries loaded at the top of the script:
```

```

## Concluding remarks
I'm working with Lee Berstrand to develop a (nearly) fully-automated method to build similar figures in a tool called BackBLAST: https://github.com/LeeBergstrand/BackBLAST_Reciprocal_BLAST. Once finished, this tool should replace the need for the instructions I've left above.


