# ABOUT data acquisition
Copyright Jackson M. Tsuji, Neufeld Research Group, 2019  
Part of the larger *IISD-ELA Chlorobia cyc2 project*.

All code here is to be run in a Bash terminal of a unix system (e.g., Ubuntu 16 or 18) unless indicated otherwise.

## Define where you downloaded the Github repo:
```bash
github_repo_location="/Analysis/jmtsuji/chlorobia-cyc2-genomics"
```

## Software prerequisites
- miniconda (miniconda3 preferred)

## Conda environment with all needed dependencies:
```bash
conda create -y -n sra_download -c bioconda sra-tools bbmap
```
Use this environment via `conda activate sra_download` (as shown below).


## Part A: Unassembled lake metagenomes
Sequencing was performed on an Illumina HiSeq 2500 with 2x200 bp reads (see Methods). Raw reads for the project can be downloaded from NCBI. The `sra_toolkit` needs to be installed, along with the BBTools suite for reformating the data. See the conda code above for this.

To download the reads (~50 GB!):
```bash
# Activate the conda environment
conda activate sra_download

# Load the info table
download_dir="${github_repo_location}/Data_analysis_pipeline/01_data_acquisition/lake_metagenomes"
lake_metagenome_info_filepath="${download_dir}/../srr_accessions_lake_metagenomes.tsv"

mkdir -p ${download_dir}
cd ${download_dir}

# Get the info for each sample
sample_IDs=($(cut -d $'\t' -f 4 ${lake_metagenome_info_filepath} | tail -n +2 ))
srr_accessions=($(cut -d $'\t' -f 5 ${lake_metagenome_info_filepath} | tail -n +2 ))

# Download the data
for i in $(seq 1 ${#srr_accessions[@]}); do
    # Set counter to zero-ordered
    j=$((${i}-1))

    # Get variables
    sample_ID=${sample_IDs[${j}]}
    accession=${srr_accessions[${j}]}

    echo "[ $(date -u) ]: Downloading '${sample_ID}' from '${accession}'"
    fastq-dump --stdout --split-spot ${accession} 2>${accession}.log | reformat.sh in=stdin.fastq int=t verifypaired=t out1="${sample_ID}_R1.fastq.gz" out2="${sample_ID}_R2.fastq.gz" 2>>${accession}.log
done

# WARNING! This might not work if the R1 and R2 reads are different lengths, for other datasets. The dataset used in this paper should be okay.
```

## Part B: Unassembled enrichment culture metagenomes
Sequencing was performed on an Illumina HiSeq 2500 with 2x125 bp reads (see Methods). Raw reads for the project can be downloaded from NCBI as above.

To download the reads:
```bash
# Activate the conda environment
conda activate sra_download

# Load the info table
download_dir="${github_repo_location}/Data_analysis_pipeline/01_data_acquisition/enrichment_metagenomes"
enrichment_info_filepath="${download_dir}/../srr_accessions_enrichment_cultures.tsv"

mkdir -p ${download_dir}
cd ${download_dir}

# Get the info for each sample
sample_IDs=($(cut -d $'\t' -f 6 ${enrichment_info_filepath} | tail -n +2 ))
srr_accessions=($(cut -d $'\t' -f 7 ${enrichment_info_filepath} | tail -n +2 ))

# Download the data
for i in $(seq 1 ${#srr_accessions[@]}); do
    # Set counter to zero-ordered
    j=$((${i}-1))

    # Get variables
    sample_ID=${sample_IDs[${j}]}
    accession=${srr_accessions[${j}]}

    echo "[ $(date -u) ]: Downloading '${sample_ID}' from '${accession}'"
    fastq-dump --stdout --split-spot ${accession} 2>${accession}.log | reformat.sh in=stdin.fastq int=t verifypaired=t out1="${sample_ID}_R1.fastq.gz" out2="${sample_ID}_R2.fastq.gz" 2>>${accession}.log
done

# WARNING! This might not work if the R1 and R2 reads are different lengths, for other datasets. The dataset used in this paper should be okay.
```
After running Parts A and B, you will be ready to run the assembly and binning pipeline. Otherwise, if you want direct access to the genome bins, see below.

## Part C: *Chlorobia* genome bins on NCBI
**NOT YET FINISHED -- enrichment culture genomes are still being annotated by NCBI**

This code can also be found as part of the analysis workflow in `06_comparative_genomics/README.md`.  

**NOTE**: You can also download the *Chlorobia* genome bins from this repo's corresponding Zenodo data repository in Part D below. Those version of the genomes are the same as what was used in the analyses for the paper. The versions on NCBI have gone through NCBI's annotation pipeline instead. There might be minor differences because of this. In general, I think it is fine to use the NCBI versions of the genomes unless you really want to exactly replicate what I did.

To download, run:
```
download_dir="${github_repo_location}/Data_analysis_pipeline/01_data_acquisition/ELA_Chlorobia_genomes"
info_filepath="${download_dir}/ELA_Chlorobia_genomes_links.tsv"
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

## Part D: *Chlorobia* genome bins and uncurated genome bins from the sequencing project
A number of additional genome bins were generated for this project that have not been manually curated. You'll see some stats on these in Figure 3 of the manuscript, for example. Although we did not deposit these on NCBI (we plan a future paper that will properly analyze and deposit the bins), they can be downloaded from a Zenodo repo.

Similar download code is also in `05_bin_analysis/README.md`.
```bash
destination_dir="${github_repo_location}/Data_analysis_pipeline/01_data_acquisition/zenodo_genome_bins"
zenodo_url="https://zenodo.org/record/3228469/files/dereplicated_genomes.tar.gz"

mkdir -p ${destination_dir}
cd ${destination_dir}

# Download
wget -O - ${zenodo_url} | tar -xzf -
# Downloads `dereplicated_genomes.tar.gz` and unpacks into `dereplicated_genomes`

# Contains two folders: curated (with *Chlorobia* bins) and uncurated (the other bins)
```

## Part E: assembled contig download
Most people do not need these, but if you are interested, you can download the entire assembled contig set of the single assemblies or co-assemblies from NCBI as nucleotide files:
```bash
download_dir="${github_repo_location}/Data_analysis_pipeline/01_data_acquisition/contigs"
info_filepath="${download_dir}/../lake_metagenome_contigs_info.tsv"
logfile="${download_dir}/lake_metagenome_contigs_info_download.log"

mkdir -p ${download_dir}
cd ${download_dir}

# Load the URLs from the provided guide file
metagenome_names=($(cut -d $'\t' -f 5 ${info_filepath} | tail -n +2))
metagenome_urls_base=($(cut -d $'\t' -f 7 ${info_filepath} | tail -n +2))

# How to work with the base URL:
# E.g., for 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/020/505/GCF_000020505.1_ASM2050v1'
# Genome is at: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/020/505/GCF_000020505.1_ASM2050v1/GCF_000020505.1_ASM2050v1_genomic.fna.gz
# No protein predictions available for contigs

echo "[ $(date -u) ]: Downloading ${#metagenome_names[@]} assembled contig sets" | tee ${logfile}

# Then download
for i in $(seq 1 ${#metagenome_names[@]}); do
# Set counter to zero-ordered
j=$((${i}-1))

metagenome_name=${metagenome_names[${j}]}
metagenome_url_base=${metagenome_urls_base[${j}]}
metagenome_ID=${metagenome_url_base##*/}

echo "[ $(date -u) ]: Downloading '${metagenome_name}' from '${metagenome_url_base}'" | tee -a ${logfile}
wget -nv -O - ${metagenome_url_base}/${metagenome_ID}_genomic.fna.gz > ${metagenome_name}.fna.gz

done
```


