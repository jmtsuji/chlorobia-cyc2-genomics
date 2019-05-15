# ABOUT comparative genomics of *Chlorobia* genome bins
Copyright Jackson M. Tsuji, Neufeld Research Group, 2019
Part of the larger *IISD-ELA Chlorobia cyc2 project*.

All code here is to be run in a Bash terminal of a unix system (e.g., Ubuntu 16 or 18) unless indicated otherwise.

Most of the analyses done here generate output files used to generate publication figures. The figure that each output data file is associated with is indicated.

## Define where you downloaded the Github repo:
```
github_repo_location="/Analysis/jmtsuji/chlorobia-cyc2-genomics"
```

## Software prerequisites
- miniconda (miniconda3 preferred)

## Conda environment with all needed dependencies:
```
conda create -y -n genome_comparison -c bioconda iqtree=1.6.10 gblocks=0.91b clustalo=1.2.3 seqtk=1.3 prodigal=2.6.3

# Install part of the basic sequence analysis suite
conda activate genome_comparison
cd /tmp
wget -nv -O - https://github.com/jmtsuji/basic-sequence-analysis/archive/v1.2.0.tar.gz | tar -xzf -
cp basic-sequence-analysis-1.2.0/basic-sequence-analysis-version basic-sequence-analysis-1.2.0/text_find_and_replace.sh ${CONDA_PREFIX}/bin
rm -rf basic-sequence-analysis-1.2.0
conda deactivate
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

# Then change extension to .fna
find . -name "*.fa" | xargs -I {} basename {} ".fa" | xargs -I {} mv {}.fa {}.fna

# Then predict amino acids via Prodigal
fna_files=($(find . -name "*.fna" | sort -h))
for fna_file in ${fna_files[@]}; do

echo "[ $(date -u) ]: Predicting genes for '${fna_file%.fna}'"
prodigal -i ${fna_file} -a ${fna_file%.fna}.faa -c -f gff -o ${fna_file%.fna}.gff 2>${fna_file%.fna}.progidal.log

done

# Gzip to match NCBI files
gzip *.fna *.faa
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
wget -nv -O - ${genome_url_base}/${genome_ID}_genomic.fna.gz > ${genome_name}.fna.gz
wget -nv -O - ${genome_url_base}/${genome_ID}_protein.faa.gz > ${genome_name}.faa.gz
wget -nv -O - ${genome_url_base}/${genome_ID}_genomic.gff.gz > ${genome_name}.gff.gz

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

Install GToTree and a helper script
```
work_dir="${github_repo_location}/Data_analysis_pipeline/06_comparative_genomics/03_chlorobia_phylogeny"
GToTree_version="1.1.10"
basic_sequence_analysis_version="1.2.2"

# Install GToTree
mkdir -p ${work_dir}/Installation_files
cd ${work_dir}/Installation_files
wget -nv -O - https://github.com/AstrobioMike/GToTree/archive/v${GToTree_version}.tar.gz > GToTree-v${GToTree_version}.tar.gz
tar -xzf GToTree-v${GToTree_version}.tar.gz
rm GToTree-v${GToTree_version}.tar.gz
cd GToTree-${GToTree_version}/
./conda-setup.sh

# Install wrapper
conda activate gtotree
cd ${work_dir}/Installation_files
wget https://github.com/jmtsuji/basic-sequence-analysis/archive/v${basic_sequence_analysis_version}.tar.gz
tar -xzf v${basic_sequence_analysis_version}.tar.gz
rm v${basic_sequence_analysis_version}.tar.gz
cd basic-sequence-analysis-${basic_sequence_analysis_version}
cp phylogeny_builder_whole_genome.sh basic-sequence-analysis-version ${CONDA_PREFIX}/bin
cd ..
rm -rf basic-sequence-analysis-${basic_sequence_analysis_version}
conda deactivate
```
Can now use this conda environment by running `conda activate gtotree`

Then make the genome phylogeny, using the `gtotree` conda env
```
work_dir="${github_repo_location}/Data_analysis_pipeline/06_comparative_genomics/03_chlorobia_phylogeny"
source_dir="${work_dir}/../01_Chlorobia_genomes"
genome_dir="${work_dir}/inputs"
output_dir="${work_dir}/phylogeny"
log_name="${work_dir}/phylogeny_builder_whole_genome.log"
threads=24
bootstrap_replicates=1000
bootstrap_type="normal"
phylogenetic_model="Universal_Hug_et_al.hmm" # run `gtt-hmms` to see all available models

mkdir -p ${genome_dir}
cd ${work_dir}

# Put the amino acid files into a common folder via hard links
find ${source_dir} -name "*.faa.gz" | xargs -I {} ln {} ${genome_dir}

# Run the whole genome phylogeny script
phylogeny_builder_whole_genome.sh -@ ${threads} -b ${bootstrap_replicates} \
    -B ${bootstrap_type} -p ${phylogenetic_model} \
    ${genome_dir} ${output_dir} 2>&1 | tee ${log_name}

```
The output file `*.treefile` can now be used in generating Figure 2 and Supplementary Figure S1.

## Make subsetted *cyc2* and riboprotein phylogenies for phylogenetic comparison
This is for Supplementary Figure 4. Subsetting to shared entries between the two sequence sets and then re-aligning and making trees. Procedure is basically the same as above with fewer sequences.

Working in `04_subset_phylogenies`

### **cyc2**
Subset the *Chlorobia* sequences with genomes, using the guide file in the repo
```
source_filepath="${github_repo_location}/Data_analysis_pipeline/04_HMM_development/01_downloads/final/cyc2_seqs_chlorobia.faa"
cyc2_info_filepath="${destination_dir}/cyc2_subset_info.tsv"
destination_dir="${github_repo_location}/Data_analysis_pipeline/06_comparative_genomics/04_subset_phylogenies/cyc2/input"
destination_filepath="${destination_dir}/cyc2_subset.faa"
mkdir -p ${destination_dir}
cd ${destination_dir}

# Get the sequences
cut -d $'\t' -f 1 ${cyc2_info_filepath} | tail -n +2 > ${destination_filepath%.faa}.list.tmp
seqtk subseq ${source_filepath} ${destination_filepath%.faa}.list.tmp > ${destination_filepath}
rm ${destination_filepath%.faa}.list.tmp

# Fix names to match the names of the genomes
text_find_and_replace.sh cyc2_subset_info.tsv cyc2_subset.faa cyc2_subset_renamed.faa
```

Align and clean with GBlocks
```
work_dir="${github_repo_location}/Data_analysis_pipeline/06_comparative_genomics/04_subset_phylogenies/cyc2"
alignment_dir="${work_dir}/align_and_clean"
input_filepath="${work_dir}/input/cyc2_subset_renamed.faa"
threads=4

mkdir -p ${alignment_dir}
cd ${alignment_dir}

# Align
base_name=${input_filepath%.faa}
base_name=${base_name##*/}
echo "[ $(date -u) ]: Aligning sequence '${base_name}.faa'"
clustalo --in ${input_filepath} --out ${base_name}_aligned.faa --distmat-out=${base_name}_distmat.txt --full --percent-id --threads=${threads} --verbose 2>&1 | tee ${base_name}_alignment.log

# Clean
input_filepath=cyc2_subset_renamed_aligned.faa
Gblocks ${input_filepath} -t=p -b3=40 -b4=4 -b5=h -e=_GB01 2>&1 | tee ${input_filepath}_GB01.log
```

Built the phylogenetic tree
```
phylogeny_dir="${work_dir}/phylogeny"
input_filepath="${alignment_dir}/cyc2_subset_renamed_aligned.faa_GB01"
threads=10

mkdir -p ${phylogeny_dir}
cd ${phylogeny_dir}

iqtree -s ${input_filepath} -pre cyc2_subset_phylogeny -nt ${threads} -seed 47 -b 1000 -m MFP
```

### riboprotein
Subset the *Chlorobia* genomes with *cyc2*
```
source_dir="${github_repo_location}/Data_analysis_pipeline/06_comparative_genomics/03_chlorobia_phylogeny/inputs"
destination_dir="${github_repo_location}/Data_analysis_pipeline/06_comparative_genomics/04_subset_phylogenies/riboprotein/inputs"
genomes_to_keep=(Chl_ferrooxidans_KoFox Chl_luteolum_DSM_273 Chl_phaeoferrooxidans_KB01 Chl_sp_N1 L227_2013_bin_56 L227_enrichment_S_6D L304_enrichment_S_6D)

mkdir -p ${destination_dir}
cd ${destination_dir}

# Make hard links of the genome files of interest
for genome_to_keep in ${genomes_to_keep[@]}; do

echo "[ $(date -u) ]: Linking '${genome_to_keep}'"
ln ${source_dir}/${genome_to_keep}.faa.gz .

done

```

Build the phylogeny, inside the `gtotree` conda env
```
work_dir="${github_repo_location}/Data_analysis_pipeline/06_comparative_genomics/04_subset_phylogenies/riboprotein"
genome_dir="${work_dir}/inputs"
output_dir="${work_dir}/phylogeny"
log_name="${work_dir}/phylogeny_builder_whole_genome.log"
threads=24
bootstrap_replicates=1000
bootstrap_type="normal"
phylogenetic_model="Universal_Hug_et_al.hmm" # run `gtt-hmms` to see all available models

cd ${work_dir}

# Run the whole genome phylogeny script
phylogeny_builder_whole_genome.sh -@ ${threads} -b ${bootstrap_replicates} \
    -B ${bootstrap_type} -p ${phylogenetic_model} \
    ${genome_dir} ${output_dir} 2>&1 | tee ${log_name}
```

Now, the tree files output by these two analyses can be cross-compared -- see the Supplementary Figure S4 folder.

## *Chlorobia* gene pathway analysis
Used for the Figure 2 heatmap
# TODO

## ANI calculation for *Chlorobia* genomes
Used for Supplementary Figure S1
# TODO

### Done!
This is the end of the main heavy-lifting data processing work for this paper. Figures were generated based off this dataset using code found in each figure folder.

