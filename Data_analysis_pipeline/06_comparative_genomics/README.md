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
conda create -y -n genome_comparison -c bioconda iqtree=1.6.10 gblocks=0.91b clustalo=1.2.3 seqtk=1.3
```
Use this environment via `conda activate genome_comparison` (as shown below).


## Download *Chlorobia* genomes
Download the *Chlorobia* genomes from this study from NCBI
```
download_dir="${github_repo_location}/Data_analysis_pipeline/06_comparative_genomics/01_Chlorobia_genomes/this_study"
info_filepath="${download_dir}/ELA_Chlorobia_links.tsv"
mkdir -p ${download_dir}
cd ${download_dir}

# TODO
# Because the NCBI accession is not yet up, download from Zenodo instead
zenodo_url=
wget -O - ${zenodo_url} | tar -xvzf -
# Makes ____
# Then keep just the Chlorobia genomes

```
The alternative would be to use the genome bins you generated yourself, if you ran the whole pipeline up to this point, at `03_bin_curation/03_contig_ordering/ordered_genomes/final`.  
If you'd like to do that, then copy those genomes into the above `download_dir` in place of the ones from NCBI.

Download the type strain *Chlorobia* genomes from NCBI
```
download_dir="${github_repo_location}/Data_analysis_pipeline/06_comparative_genomics/01_Chlorobia_genomes/reference"
info_filepath="${download_dir}/reference_Chlorobia_links.tsv"
logfile="${download_dir}/download.log"

mkdir -p ${download_dir}
cd ${download_dir}

# Load the URLs from the provided guide file
genome_names=($(cut -d $'\t' -f 1 ${info_filepath} | tail -n +2))
genome_urls_base=($(cut -d $'\t' -f 3 ${info_filepath} | tail -n +2))

# How to work with the base URL:
# E.g., for 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/020/505/GCF_000020505.1_ASM2050v1'
# Genome is at: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/020/505/GCF_000020505.1_ASM2050v1/GCF_000020505.1_ASM2050v1_genomic.fna.gz
# Proteins are at: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/020/505/GCF_000020505.1_ASM2050v1/GCF_000020505.1_ASM2050v1_protein.faa.gz
# GFF is at: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/020/505/GCF_000020505.1_ASM2050v1/GCF_000020505.1_ASM2050v1_genomic.gff.gz

echo "[ $(date -u) ]: Downloading ${#genome_names[@]} genomes" | tee ${logfile}

# Then download
for i in $(seq 1 ${#genome_names[@]}); do
# Set counter to zero-ordered
j=$((${i}-1))

genome_name=${genome_names[${j}]}
genome_url_base=${genome_urls_base[${j}]}
genome_ID=${genome_url_base##*/}

echo "[ $(date -u) ]: Downloading '${genome_name}' from '${genome_url_base}'" | tee -a ${logfile}
wget -nv -O - ${genome_url_base}/${genome_ID}_genomic.fna.gz > ${genome_name}.fna
wget -nv -O - ${genome_url_base}/${genome_ID}_protein.faa.gz > ${genome_name}.faa
wget -nv -O - ${genome_url_base}/${genome_ID}_genomic.gff.gz > ${genome_name}.gff

done
```

# TODO - move this whole section to Figure 1 instead!!!
## Gene neighbourhood of *Chlorobia* cyc2
Used for Figure 1, panel B

Note that the accessions of *cyc2* genes within each applicable genome are summarized in `02_cyc2_gene_neighbourhood/Chlorobia_cyc2_genome_info.tsv`
# TODO - add accessions for the ELA genome bins!!!

Pull out the neighbouring genes around the *cyc2* gene for each genome using a custom R script
# TODO - move this whole section to Figure 1 instead!!!


## Phylogeny of *cyc2* genes
Used for Figure 1, panel C

Aligned collection of *cyc2* is already available from `04_HMM_development`. Copy it into this folder:
```
work_dir="${github_repo_location}/Data_analysis_pipeline/06_comparative_genomics/02_cyc2_phylogeny"
Gblocks_dir="${work_dir}/01_Gblocks"
source_filepath="${github_repo_location}/Data_analysis_pipeline/04_HMM_development/02_alignment/cyc2_seqs_all_aligned.faa"
mkdir -p ${Gblocks_dir}
cd ${Gblocks_dir}

cp ${source_filepath} ${Gblocks_dir}
```

Masked the alignment using Gblocks `v0.91b`:
```
input_filepath="cyc2_seqs_all_aligned.faa"
Gblocks ${input_filepath} -t=p -b3=40 -b4=4 -b5=h -e=_GB01 2>&1 | tee ${input_filepath}_GB01.log
```

Built the phylogenetic tree via IQ-TREE `v1.6.10`:
```
phylogeny_dir="${work_dir}/02_phylogeny"
input_filepath="${Gblocks_dir}/cyc2_seqs_all_aligned.faa_GB01"
threads=10

mkdir -p ${phylogeny_dir}
cd ${phylogeny_dir}

iqtree -s ${input_filepath} -pre cyc2_phylogeny -nt ${threads} -seed 47 -b 1000 -m MFP
```
The output file `cyc2_phylogeny.contree` contains the consensus tree data for the phylogenetic tree.  

See Figure 1 folder for the code of how this was plotted.


## Phylogenetic placement of *Chlorobia* genome bins
Used for Figure 2

```
work_dir="${github_repo_location}/Data_analysis_pipeline/06_comparative_genomics"
source_dir="${github_repo_location}/Data_analysis_pipeline/06_comparative_genomics/01_Chlorobia_genomes"
genome_dir="${work_dir}/03_chlorobia_phylogeny/inputs"
output_dir="${work_dir}/03_chlorobia_phylogeny/phylogeny"
log_name="${work_dir}/phylogeny_builder_whole_genome.log"
threads=24
bootstrap_replicates=1000
bootstrap_type="normal"
phylogenetic_model="Universal_Hug_et_al.hmm" # run `gtt-hmms` to see all available models

# Put the amino acid files into a common folder via hard links
find ${source_dir} -iname "*.faa.gz" | xargs -I {} ln {} ${genome_dir}

# Run the script
mkdir -p ${work_dir}
cd ${work_dir}
phylogeny_builder_whole_genome.sh -@ ${threads} -b ${bootstrap_replicates} \
    -B ${bootstrap_type} -p ${phylogenetic_model} \
    ${genome_dir} ${output_dir} 2>&1 | tee ${log_name}

```
The output file `*.treefile` can now be used in generating Figure 2 and Supplementary Figure S1.



OLD

### Identify shared ribosomal proteins
This was done by manually looking through annotation files of each genome. 
It was determined that nine key ribosomal protein genes are shared by all genomes (see manuscript). 
The IDs of these genes are summarized in `${github_repo_location}/Data_analysis_pipeline/06_comparative_genomics/03_chlorobia_phylogeny/01_riboproteins/shared_riboproteins.tsv`

### Extract and align the shared proteins
```
# User variables
genome_dir="${github_repo_location}/Data_analysis_pipeline/06_comparative_genomics/01_Chlorobia_genomes"
output_dir="${github_repo_location}/Data_analysis_pipeline/06_comparative_genomics/03_chlorobia_phylogeny"

# Load table data
genome_names=($(cut -d $'\t' -f 1 ${guide_filepath} | tail -n +2))
gene_names=($(cut -d $'\t' -f 2 ${guide_filepath} | tail -n +2))
gene_accessions=($(cut -d $'\t' -f 3 ${guide_filepath} | tail -n +2))

# Output directory cannot already exist
if [ -d ${output_dir} ]; then
echo "[ $(date -u) ]: Cannot continue because '${output_dir}' already exists. Please delete the output directory and then try again."
else

for i in $(seq 1 ${#gene_name[@]}); do
# Set counter to zero-order
j=$((${i}-1))

# Get variables
gene_name=${gene_names[${j}]}
genome_name=${genome_names[${j}]}
gene_accession=${gene_accessions[${j}]}

done
fi

```

### Generate ribosomal protein phylogeny


## *Chlorobia* gene pathway analysis
Used for Figure 2 heatmap

## Comparative phylogeny of *cyc2* and ribosomal proteins
Used for Supplementary Figure S4



### Done!
This is the end of the main heavy-lifting data processing work for this paper. Figures were generated based off this dataset using code found in each figure folder.

