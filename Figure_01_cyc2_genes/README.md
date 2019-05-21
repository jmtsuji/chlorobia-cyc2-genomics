# ABOUT Figure 1 - cyc2 gene verification and novelty
Copyright Jackson M. Tsuji, Neufeld Research Group, 2019
Part of the larger *IISD-ELA Chlorobia cyc2 project*.

Each panel of this figure was generated independently. All panels were then assembled manually.

## 1. Data collection
- Novel _cyc2_ from this project: predicted open reading frames (ORFs) from all assembled contigs from this project were scanned using MetAnnotate for _Chlorobia_ cyc2 sequences using a custom HMM developed for this project (see Methods). Three of the novel _cyc2_ were contained in high-quality dereplicated _Chlorobia_ genome bins (see ABOUT files for assembly and binning). The last was contained in a medium-quality dereplicated genome bin. One additional potential _cyc2_ sequence belonging to _Chlorobia_ was found but was excluded because of lacking the c5 family cytochrome gene neighbour (it was on a short contig).
- Reference _cyc2_ sequences: prediced ORFs from all reference organisms containing _cyc2_ homologues described by Hu and colleagues (2017) were downloaded from NCBI. A BLASTP search using the _Mariprofundus ferroxydans_ PV-1 _cyc2_ homologue was used to identify _cyc2_ homologues in the other genomes.
- Gene neighbourhood information: GFF3 files were obtained for all genome bins (or contig collections) from this project. GFF3 files were also downloaded from NCBI for the reference _Chlorobia_ strains known to possess _cyc2_ homologues.

Code to download reference organism genomes, ORF predictions, and GFF3 files from NCBI:
```

```

## 2. Panel 1: conceptual model of the Cyc2 protein
No special code here. As briefly described in the paper, the _cyc2_ predicted amino acid sequence for _Chlorobium phaeoferrooxidans_ KB01 was uploaded to the Phyre2 webserver to predict the folding structure based on homology to known PDB entries. The resulting PDB file was visualized in Chimera2.

## 3. Panel 2: _Chlorobia_ _cyc2_ gene neighbour analysis
The GFF3 files of all input bins or contig collections are available in the `04_extra_files` folder in this repo.

```
# TODO - move this whole section to Figure 1 instead!!!
## Gene neighbourhood of *Chlorobia* cyc2
Used for Figure 1, panel B

Note that the accessions of *cyc2* genes within each applicable genome are summarized in `02_cyc2_gene_neighbourhood/Chlorobia_cyc2_genome_info.tsv`
# TODO - add accessions for the ELA genome bins!!!

Pull out the neighbouring genes around the *cyc2* gene for each genome using a custom R script
# TODO - move this whole section to Figure 1 instead!!!
```



The following R code was used to generate the base visualization, taking the GFF3 files as input. You'll need to adjust the working directory to match whatever folder you work in. You'll also need to install all libraries loaded at the top of the script:
```

```

After this, the figure was a bit messy and incomplete. Manual annotation work was used to clean up the figure and add additional colours and details. Some of this could be coded long-term (I would welcome this!), but it exceeded the scope of what was possible for this project.

Specifically, the following changes to the base figure were done:
- Names of the organisms (left side) were simplified, and the font was made larger.
- Labels for all gene names were (painstakingly) added by browsing through the GFF files in Excel, finding the genes correponding to each entry, searching the name of the predicted protein in InterPro, NCBI, and the literature, and choosing an appropriate gene abbreviation.
- Closest homologues of all genes were (painstakingly) checked via BLASTP. The predicted ORFs of the flanking ~50 genes on either side of the predicted _cyc2_ gene were obtained from the `.faa` files mentioned in the Data Collection section above and were searched via BLASTP. See methodology for classifying gene types in the Figure 2 caption of the paper.
- Unique colours were added for genes with the same predicted product between different organisms and for genes that had a closest homologe classified to a non-_Chlorobia_ organism.
- A colour legend was added at the bottom of the figure.

## 4. Panel 3: _cyc2_ phylogeny and multiple sequence alignemnt
### Part A -- alignment and phylogeny
The following code was run:
```

```

### Part B -- producing the figure
Producing the figure was almost entired automated. The following input files are needed:
- Non-masked _cyc2_ sequence alignment generated above
- _cyc2_ gene phylogeny generated above

The following R code was used to generate the figure. You'll need to adjust the working directory to match whatever folder you work in. You'll also need to install all libraries loaded at the top of the script:
```

```
Any post-production edits to this figure were minimal (e.g., making fonts bold, simplying label names, highlighting sequences with colours).


