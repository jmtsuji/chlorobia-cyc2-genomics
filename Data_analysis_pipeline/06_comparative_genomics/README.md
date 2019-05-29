# ABOUT comparative genomics of *Chlorobia* genome bins
Copyright Jackson M. Tsuji, Neufeld Research Group, 2019  
Part of the larger *IISD-ELA Chlorobia cyc2 project*.

All code here is to be run in a Bash terminal of a unix system (e.g., Ubuntu 16 or 18) unless indicated otherwise.

Most of the analyses done here generate output files used to generate publication figures. The figure that each output data file is associated with is indicated.

## Define where you downloaded the Github repo:
```bash
github_repo_location="/Analysis/jmtsuji/chlorobia-cyc2-genomics"
```

## Software prerequisites
- miniconda (miniconda3 preferred)

## Conda environment with needed dependencies for the majority of the README:
```bash
conda create -y -n genome_comparison -c bioconda iqtree=1.6.10 gblocks=0.91b clustalo=1.2.3 seqtk=1.3 prodigal=2.6.3 fastani=1.1

# Install part of the basic sequence analysis suite
conda activate genome_comparison
cd /tmp
wget -nv -O - https://github.com/jmtsuji/basic-sequence-analysis/archive/v1.2.0.tar.gz | tar -xzf -
cp basic-sequence-analysis-1.2.0/basic-sequence-analysis-version basic-sequence-analysis-1.2.0/text_find_and_replace.sh ${CONDA_PREFIX}/bin
rm -rf basic-sequence-analysis-1.2.0
conda deactivate
```
Use this environment via `conda activate genome_comparison` (as shown below).  
This is used for all analyses except for the BackBLAST core module and GToTree, which have their own conda envs (as shown below)

## Download *Chlorobia* genomes

### Download the *Chlorobia* genomes from this study
You have several options:

**Option 1** (recommended): Download from the Zenodo repo. This is the method I've tested and know works.
```bash
download_dir="${github_repo_location}/Data_analysis_pipeline/06_comparative_genomics/01_Chlorobia_genomes/this_study"
zenodo_url="https://zenodo.org/record/3228469/files/dereplicated_genomes.tar.gz"

mkdir -p ${download_dir}
cd ${download_dir}

# Download
wget -O - ${zenodo_url} | tar -xzf -
# Downloads `dereplicated_genomes.tar.gz` and unpacks into `dereplicated_genomes`

# Change extension to .fna to look like NCBI
find dereplicated_genomes/curated -type f -name "*.fa" | xargs -I {} basename {} ".fa" | xargs -I {} mv {}.fa {}.fna

# Move the curated bins all into the main folder, then delete the non-curated ones
find dereplicated_genomes/curated -type f -name "*.fa" | xargs -I {} mv {} .
rm -r dereplicated_genomes

# Then predict amino acids via Prodigal
fna_files=($(find . -name "*.fna" | sort -h))
for fna_file in ${fna_files[@]}; do

    echo "[ $(date -u) ]: Predicting genes for '${fna_file%.fna}'"
    prodigal -i ${fna_file} -a ${fna_file%.fna}.faa -c -f gff -o ${fna_file%.fna}.gff 2>${fna_file%.fna}.progidal.log

done

# Gzip to match NCBI files
gzip *.fna *.faa *.gff
```

**Option 2**: download from NCBI. This is a nice option, but the genomes were not yet available on NCBI at the time of writing, so I am not sure how NCBI's annotations will compare to those run on our lab's server. There might be a few differences. Also, the code has not been tested to work with these NCBI genome files end-to-end, so there could be bugs.  
**NOT YET FINISHED -- enrichment culture genomes are still being annotated by NCBI**
```bash
download_dir="${github_repo_location}/Data_analysis_pipeline/06_comparative_genomics/01_Chlorobia_genomes/this_study"
info_filepath="${download_dir}/ELA_Chlorobia_links.tsv"
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

**Option 3**: use the genome bins you generated yourself, if you ran the whole pipeline up to this point, at `03_bin_curation/03_contig_ordering/ordered_genomes/final`.  
If you'd like to do that, then copy those genomes into the above `download_dir` in place of the ones from NCBI.

### Download the type strain (reference) *Chlorobia* genomes from NCBI
```bash
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

## Phylogeny of *cyc2* genes
Used for Figure 1, panel C

Aligned collection of *cyc2* is already available from `04_HMM_development`. Copy it into this folder:
```bash
work_dir="${github_repo_location}/Data_analysis_pipeline/06_comparative_genomics/02_cyc2_phylogeny"
Gblocks_dir="${work_dir}/01_Gblocks"
source_filepath="${github_repo_location}/Data_analysis_pipeline/04_HMM_development/02_alignment/cyc2_seqs_all_aligned.faa"
mkdir -p ${Gblocks_dir}
cd ${Gblocks_dir}

cp ${source_filepath} ${Gblocks_dir}
```

Masked the alignment using Gblocks `v0.91b`:
```bash
input_filepath="cyc2_seqs_all_aligned.faa"
Gblocks ${input_filepath} -t=p -b3=40 -b4=4 -b5=h -e=_GB01 2>&1 | tee ${input_filepath}_GB01.log
```

Built the phylogenetic tree via IQ-TREE `v1.6.10`:
```bash
phylogeny_dir="${work_dir}/02_phylogeny"
input_filepath="${Gblocks_dir}/cyc2_seqs_all_aligned.faa_GB01"
threads=10

mkdir -p ${phylogeny_dir}
cd ${phylogeny_dir}

iqtree -s ${input_filepath} -pre cyc2_phylogeny -nt ${threads} -seed 47 -b 1000 -m MFP
```
The output file `cyc2_phylogeny.treefile` contains the phylogenetic tree info.  

See Figure 1 folder for the code of how this was plotted.


## Phylogenetic placement of *Chlorobia* genome bins
Used for Figure 2

Install GToTree and a helper script
```bash
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
```bash
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
This is for Supplementary Figure 3. Subsetting to shared entries between the two sequence sets and then re-aligning and making trees. Procedure is basically the same as above with fewer sequences.

Working in `04_subset_phylogenies`

### **cyc2**
Subset the *Chlorobia* sequences with genomes, using the guide file in the repo
```bash
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
```bash
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
```bash
phylogeny_dir="${work_dir}/phylogeny"
input_filepath="${alignment_dir}/cyc2_subset_renamed_aligned.faa_GB01"
threads=10

mkdir -p ${phylogeny_dir}
cd ${phylogeny_dir}

iqtree -s ${input_filepath} -pre cyc2_subset_phylogeny -nt ${threads} -seed 47 -b 1000 -m MFP
```

### riboprotein
Subset the *Chlorobia* genomes with *cyc2*
```bash
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
```bash
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

Now, the tree files output by these two analyses can be cross-compared -- see the Supplementary Figure 3 folder.

## *Chlorobia* gene pathway analysis
Used for the Figure 2 heatmap

Installed BackBLAST development commit `cdaddc4`
```
work_dir="${github_repo_location}/Data_analysis_pipeline/06_comparative_genomics/05_pathway_analysis"

mkdir -p ${work_dir}
cd ${work_dir}

# First installed the dependencies via conda
conda create -y -n backblast_core -c bioconda python=2.7 biopython=1.72 blast=2.6.0

# Then got the script from the repo
git clone https://github.com/LeeBergstrand/BackBLAST_Reciprocal_BLAST.git
cd BackBLAST_Reciprocal_BLAST
git checkout cdaddc4
chmod 755 BackBLAST.py
cd ..
# Script can now be run locally via BackBLAST_Reciprocal_BLAST/BackBLAST.py
```
Run via `conda activate backblast_core`

Collected input .faa files
```
work_dir="${github_repo_location}/Data_analysis_pipeline/06_comparative_genomics/05_pathway_analysis"
source_dir="${github_repo_location}/Data_analysis_pipeline/06_comparative_genomics/01_Chlorobia_genomes"
genome_dir="${work_dir}/01_genomes"

mkdir -p ${genome_dir}
cd ${work_dir}

# Put the amino acid files into a common folder (must unzip)
faa_files=($(find ${source_dir} -name "*.faa.gz" | sort -h))

for faa_file in ${faa_files[@]}; do

	faa_basename=${faa_file##*/}
	faa_basename=${faa_basename%.faa.gz}

	echo "[ $(date -u) ]: Copying ${faa_file##*/}"
	cat ${faa_file} | gunzip -c > ${genome_dir}/${faa_basename}.faa

done
```

The manually determined query sequences are already available in this repo in the folder `02_queries`:
- `Chl_ferrooxidans_KoFox_gene_targets`: query genes selected from `Chl_ferrooxidans_KoFox.faa`
- `Chl_clathratiforme_BU_1_gene_targets.faa`: query genes selected from `Chl_clathratiforme_BU_1_gene_targets.faa`

Ran the core BackBLAST module against all inputs
```bash
work_dir="${github_repo_location}/Data_analysis_pipeline/06_comparative_genomics/05_pathway_analysis"
genome_dir="${work_dir}/01_genomes"
query_dir="${work_dir}/02_queries"
output_dir="${work_dir}/03_backblast"
log_name="${work_dir}/backblast.log"
e_value=1e-40 # maximum e-value allowable
identity=20 # minimum percent identity of hits

mkdir -p ${output_dir}/logs
cd ${work_dir}
printf "" > ${log_name}

# Find relevant files
query_files=($(find ${query_dir} -iname "*.faa" | sort -h))
ORF_files=($(find ${genome_dir} -iname "*.faa" | sort -h))

# Run BackBLAST
echo "[ $(date -u) ]: Running BackBLAST core module" | tee -a ${log_name}
echo "[ $(date -u) ]: Settings: e value cutoff = '${e_value}'; percent identity cutoff = ${identity}" | tee -a ${log_name}

for query in ${query_files[@]}; do
	query_basename=${query##*/}
	query_basename=${query_basename%.faa}

	# Get the name of the corresponding whole genome .faa file (HARD-CODED method based on file naming structure)
	query_genome_basename=${query_basename%_gene_targets}
	query_genome="${genome_dir}/${query_genome_basename}.faa"

	echo "[ $(date -u) ]: Searching for '${query_basename}' among the subjects. Using '${query_genome_basename}' as the reference predicted proteome." | tee -a ${log_name}

	for subject in ${ORF_files[@]}; do
		subject_basename=${subject##*/}
		subject_basename=${subject_basename%.faa}

		output_filename="${output_dir}/${query_genome_basename}__to__${subject_basename}.csv"
		output_logname="${output_dir}/logs/${query_genome_basename}__to__${subject_basename}.log"

		echo "[ $(date -u) ]: 'BackBLAST.py -q ${query##*/} -r ${query_genome##*/} -s ${subject##*/} -e ${e_value} -i ${identity} -o ${output_filename##*/} 2>&1 | tee ${output_logname##*/}'" | tee -a ${log_name}

		${work_dir}/BackBLAST_Reciprocal_BLAST/BackBLAST.py -q ${query} -r ${query_genome} -s ${subject} -e ${e_value} \
		    -i ${identity} -o ${output_filename} 2>&1 | tee ${output_logname}

	done
done
echo "[ $(date -u) ]: Finished." | tee -a ${log_name}
rm tempQuery.faa
```

Replaced *C. ferro* and *C. clathratiforme* with one-way BLAST due to issue with BackBLAST when BLAST'ing to self
When comparing a query to its own genome, BackBLAST currently seems to omit hits to the original gene -- maybe because of how the graph is constructed. This means that special measures are needed to compare a genome to itself. Should be pretty rudimentary; all hits should return 100% to self. No need to do reciprocal blast if the table looks normal.

Ran one-way BLAST on both queries to themselves
```bash
mkdir -p ${output_dir}/unused

query_file=${query_dir}/Chl_ferrooxidans_KoFox_gene_targets.faa
ORF_file=${genome_dir}/Chl_ferrooxidans_KoFox.faa
output_filepath="${output_dir}/Chl_ferrooxidans_KoFox__to__Chl_ferrooxidans_KoFox.csv" # same name as the current file; will replace.
e_value=1e-40
identity=20 # Not actually needed here, but I kept this so I would remember to set the pident threshold later

mv ${output_filepath} ${output_dir}/unused
blastp -query ${query_file} -subject ${ORF_file} -evalue ${e_value} -outfmt "10 qseqid sseqid pident evalue qcovhsp bitscore" > ${output_filepath}

query_file=${query_dir}/Chl_clathratiforme_BU_1_gene_targets.faa
ORF_file=${genome_dir}/Chl_clathratiforme_BU_1.faa
output_filepath="${output_dir}/Chl_clathratiforme_BU_1__to__Chl_clathratiforme_BU_1.csv" # same name as the current file; will replace.
e_value=1e-40
identity=20 # Not actually needed here, but I kept this so I would remember to set the pident threshold later

mv ${output_filepath} ${output_dir}/unused
blastp -query ${query_file} -subject ${ORF_file} -evalue ${e_value} -outfmt "10 qseqid sseqid pident evalue qcovhsp bitscore" > ${output_filepath}
```
The `.csv` output files in the main folder of `03_backblast` are summarized and combined with the *Chlorobi* phylogenetic tree generated above to produce Figure 2 -- see the Figure 2 folder.

## ANI calculation for *Chlorobia* genomes
Used for Supplementary Figure 1

```bash
work_dir="${github_repo_location}/Data_analysis_pipeline/06_comparative_genomics/06_ANI"
source_dir="$
{github_repo_location}/Data_analysis_pipeline/06_comparative_genomics/01_Chlorobia_genomes"
genome_dir="${work_dir}/inputs"
output_filename="Chlorobia_FastANI_results.txt"
log_name="${work_dir}/fastani.log"
threads=2

mkdir -p ${genome_dir}
cd ${work_dir}

# To be easy to find for future reference, put the genome files into a common folder via hard links
find ${source_dir} -name "*.fna.gz" | xargs -I {} ln {} ${genome_dir}

# Made reference and query lists for FastANI; both are the same because I am doing an all-by-all comparison
find ${genome_dir} -iname "*.fna.gz" > genomes.list

# Ran FastANI
fastANI --rl genomes.list --ql genomes.list -o ${output_filename} -t ${threads} 2>&1 | tee ${log_name}
# Runs very quickly
```

The output file `Chlorobia_FastANI_results.txt` is used in Supplementary Figure 1 (see that folder for the scripting details).


## Alignment of cytochrome c5 primary sequences
The *cyc2* gene in *Chlorobia* appears to consistently be adjacent to a c5 family cytochrome in the genome. I gathered and aligned these sequences to include in Supplementary File 5.

Obtained reference sequences manually based on the info in `reference_Chlorobia_c5_protein_info.tsv`. **MANUALLY** added sequences from this study. (Didn't code this - just a quick side analysis. Looked for c5 genes immediately adjacent to the cyc2 genes, as I already knew were there based on Figure 1B.) This resulted in `c5_family_Chlorobia_unaligned.faa`

Then aligned the c5 sequences
```bash
cd "${github_repo_location}/Data_analysis_pipeline/06_comparative_genomics/07_c5_alignment"

clustalo -i c5_family_Chlorobia_unaligned.faa --full --percent-id --distmat-out=c5_family_Chlorobia_aligned_distmat.txt -o c5_family_Chlorobia_aligned.faa --threads=2 --verbose 2>&1 | tee c5_family_Chlorobia_aligned.log
```
A copy of `c5_family_Chlorobia_aligned.faa` is included in Supplementary File 5.


## Obtaining subsets of GFF files around cyc2
For genomic context analysis in Figure 1B. Obtaining subsets makes the GFFs much easier to work with than when having full-sized files.

For reference genomes (downloaded GFF)
```
work_dir="${github_repo_location}/Data_analysis_pipeline/06_comparative_genomics/08_gff_subset/reference"
reference_chlorobia_guide_file="${work_dir}/../reference_Chlorobia_cyc2_info.tsv"
logfile="${work_dir}/reference_Chlorobia_cyc2_gff_subset.log"

mkdir -p ${work_dir}
cd ${work_dir}
printf "" > ${logfile}
echo "[ $(date -u) ]: Subsetting cyc2 for $(($(cat ${ELA_chlorobia_guide_file} | wc -l)-1)) entries" | tee -a ${logfile}

for i in $(seq 1 $(($(cat ${reference_chlorobia_guide_file} | wc -l)-1))); do
# Set counter to start at 2 to skip header
j=$((${i}+1))

# Assign variables
genome_ID=$(cut -d $'\t' -f 1 ${reference_chlorobia_guide_file} | tail -n +${j} | head -n 1)
cyc2_accession=$(cut -d $'\t' -f 2 ${reference_chlorobia_guide_file} | tail -n +${j} | head -n 1)
gff_link=$(cut -d $'\t' -f 3 ${reference_chlorobia_guide_file} | tail -n +${j} | head -n 1)

echo "[ $(date -u) ]: Downloading '${genome_ID}' from '${gff_link}'" | tee -a ${logfile}

# Download gff
wget -q -O - ${gff_link} > ${genome_ID}.gff.gz

echo "[ $(date -u) ]: Searching for protein accession '${cyc2_accession}'" | tee -a ${logfile}
matching_contig=$(zgrep "protein_id=${cyc2_accession}" ${genome_ID}.gff.gz | cut -d $'\t' -f 1)
if [ ${#matching_contig[@]} != 1 ]; then
echo "[ $(date -u) ]: Ran into problems finding '${cyc2_accession}'. Exiting..." | tee -a ${logfile}
#exit 1
fi

echo "[ $(date -u) ]: Found on contig '${matching_contig}'" | tee -a ${logfile}
zgrep "^${matching_contig}" ${genome_ID}.gff.gz > ${genome_ID}_cyc2_contig.gff
echo "[ $(date -u) ]: Grabbed $(cat ${genome_ID}_cyc2_contig.gff | wc -l) total lines corresponding to that contig. Saved as '${genome_ID}_cyc2_contig.gff'" | tee -a ${logfile}

if [ $(cat ${genome_ID}_cyc2_contig.gff | wc -l) -gt 401 ]; then
echo "[ $(date -u) ]: Output file is too large (e.g., due to long contig). Truncating to +/- 200 entries from cyc2 (same output file)" | tee -a ${logfile}
zgrep "^${matching_contig}" ${genome_ID}.gff.gz | grep -A 200 -B 200 "protein_id=${cyc2_accession}" > ${genome_ID}_cyc2_contig.gff
echo "[ $(date -u) ]: Left with $(cat ${genome_ID}_cyc2_contig.gff | wc -l) total lines. (If the number is less than 401, then it could be that the cyc2 is near the end of the contig on one side" | tee -a ${logfile}
fi

# Clean up
rm ${genome_ID}.gff.gz

done
```

For ELA *Chlorobia* cyc2 (worked directly from ATLAS output - you won't be able to replicate this unless you run the whole assembly pipeline, unfortunately)
```
work_dir="${github_repo_location}/Data_analysis_pipeline/06_comparative_genomics/08_gff_subset/ELA"
atlas_dir="${github_repo_location}/Data_analysis_pipeline/02_assembly_and_binning"
ELA_chlorobia_guide_file="${work_dir}/../ELA_Chlorobia_cyc2_info.tsv"
logfile="${work_dir}/ELA_Chlorobia_cyc2_gff_subset.log"

printf "" > ${logfile}
mkdir -p ${work_dir}
cd ${work_dir}

echo "[ $(date -u) ]: Subsetting cyc2 for $(($(cat ${ELA_chlorobia_guide_file} | wc -l)-1)) entries" | tee -a ${logfile}
for i in $(seq 1 $(($(cat ${ELA_chlorobia_guide_file} | wc -l)-1))); do
# Set counter to start at 2 to skip header
j=$((${i}+1))

# Assign variables
genome_ID=$(cut -d $'\t' -f 1 ${ELA_chlorobia_guide_file} | tail -n +${j} | head -n 1)
cyc2_accession=$(cut -d $'\t' -f 2 ${ELA_chlorobia_guide_file} | tail -n +${j} | head -n 1)
gff_filename=$(cut -d $'\t' -f 3 ${ELA_chlorobia_guide_file} | tail -n +${j} | head -n 1)

echo "[ $(date -u) ]: '${genome_ID}': finding cyc2 '${cyc2_accession}' in '${gff_filename}'" | tee -a ${logfile}

# Find the gff
gff_filepath=($(find -L ${atlas_dir} -name ${gff_filename} | grep "annotation/prokka"))
if [ ${#gff_filepath[@]} != 1 ]; then
echo "[ $(date -u) ]: Ran into problems finding '${gff_filename}'. Exiting..." | tee -a ${logfile}
#exit 1
else
echo "[ $(date -u) ]: Found GFF at '${gff_filepath}'" | tee -a ${logfile}
fi

echo "[ $(date -u) ]: Searching for protein accession '${cyc2_accession}'" | tee -a ${logfile}
matching_contig=($(grep "locus_tag=${cyc2_accession}" ${gff_filepath} | cut -d $'\t' -f 1))
if [ ${#matching_contig[@]} != 1 ]; then
echo "[ $(date -u) ]: Ran into problems finding '${cyc2_accession}'. Exiting..." | tee -a ${logfile}
#exit 1
fi

echo "[ $(date -u) ]: Found on contig '${matching_contig}'" | tee -a ${logfile}
grep "^${matching_contig}" ${gff_filepath} > ${genome_ID}_cyc2_contig.gff
echo "[ $(date -u) ]: Grabbed $(cat ${genome_ID}_cyc2_contig.gff | wc -l) total lines corresponding to that contig. Saved as '${genome_ID}_cyc2_contig.gff'" | tee -a ${logfile}

if [ $(cat ${genome_ID}_cyc2_contig.gff | wc -l) -gt 401 ]; then
echo "[ $(date -u) ]: Output file is too large (e.g., due to long contig). Truncating to +/- 200 entries from cyc2 (same output file)" | tee -a ${logfile}
grep "^${matching_contig}" ${gff_filepath} | grep -A 200 -B 200 "locus_tag=${cyc2_accession}" > ${genome_ID}_cyc2_contig.gff
echo "[ $(date -u) ]: Left with $(cat ${genome_ID}_cyc2_contig.gff | wc -l) total lines. (If the number is less than 401, then it could be that the cyc2 is near the end of the contig on one side" | tee -a ${logfile}
fi

done
```
These GFF subset files are used as input for generating Figure 1, Panel B. See that folder for details.

## Done!
This is the end of the main heavy-lifting data processing work for this paper. Figures were generated based off this dataset using code found in each figure folder.

