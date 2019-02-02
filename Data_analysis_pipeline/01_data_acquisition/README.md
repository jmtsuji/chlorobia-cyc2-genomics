# ABOUT data acquisition
Copyright Jackson M. Tsuji, Neufeld Research Group, 2019
Part of the larger *IISD-ELA Chlorobia cyc2 project*.

All code here is to be run in a Bash terminal of a unix system (e.g., Ubuntu 16 or 18) unless indicated otherwise.

## Define where you downloaded the Zenodo repo:
```
zenodo_repo_location="/Analysis/jmtsuji/Chlorobia_cyc2_code"
```

## Software prerequisites
- miniconda (miniconda3 preferred)


## Conda environment with all needed dependencies:
```
conda create -y -n sra_download -c bioconda sra-tools bbmap
```
Use this environment via `conda activate sra_download` (as shown below).


## Part A: Lake metagenomes
### 1. Data access
Sequencing was performed on an Illumina HiSeq 2500 with 2x200 bp reads (see Methods). Raw reads for the project can be downloaded from NCBI. The `sra_toolkit` needs to be installed, along with the BBTools suite for reformating the data. See the conda code above for this.

#### Download the reads (~50 GB!)
```
# Activate the conda environment
conda activate sra_download

# Load the info table
download_info_dir=${zenodo_repo_location}/Data_analysis_pipeline/01_data_acquisition
cd ${download_info_dir}
lake_metagenome_info_filename="srr_accessions_lake_metagenomes.tsv"

# Get the info for each sample
sample_IDs=($(cut -d '\t' -f 4 ${lake_metagenome_info_filename} | tail -n +2 ))
srr_accessions=($(cut -d '\t' -f 5 ${lake_metagenome_info_filename} | tail -n +2 ))

# Download the data
mkdir -p ${download_info_dir}/lake_metagenomes
cd ${download_info_dir}/lake_metagenomes

for i in $(seq 1 ${#srr_accessions[@]}); do
# Set counter to zero-ordered
j=$((${i}-1))

# Get variables
sample_ID=${sample_IDs[${j}]}
accession=${srr_accessions[${j}]}

echo "[ $(date -u) ]: Downloading '${sample_ID}' from '${accession}'"
fastq-dump --stdout --split-spot ${accession} 2>${accession}.log | reformat.sh in=stdin.fastq int=t verifypaired=t out1="${sample_ID}_R1.fastq.gz" out2="${sample_ID}_R2.fastq.gz" 2>>${accession}.log
done

# WARNING! This might not work if the R1 and R2 reads are different lengths... but I know that my data is okay.
```


## Part B: Enrichment cultures
### 1. Data access
Sequencing was performed on an Illumina HiSeq 2500 with 2x125 bp reads (see Methods). Raw reads for the project can be downloaded from NCBI as above.

#### Download the reads
```
# Activate the conda environment
conda activate sra_download

# Load the info table
download_info_dir=${zenodo_repo_location}/Data_analysis_pipeline/01_data_acquisition
cd ${download_info_dir}
enrichment_info_filename="srr_accessions_enrichment_cultures.tsv"

# Get the info for each sample
sample_IDs=($(cut -d '\t' -f 6 ${enrichment_info_filename} | tail -n +2 ))
srr_accessions=($(cut -d '\t' -f 7 ${enrichment_info_filename} | tail -n +2 ))

# Download the data
mkdir -p ${download_info_dir}/enrichment_metagenomes
cd ${download_info_dir}/enrichment_metagenomes

for i in $(seq 1 ${#srr_accessions[@]}); do
# Set counter to zero-ordered
j=$((${i}-1))

# Get variables
sample_ID=${sample_IDs[${j}]}
accession=${srr_accessions[${j}]}

echo "[ $(date -u) ]: Downloading '${sample_ID}' from '${accession}'"
fastq-dump --stdout --split-spot ${accession} 2>${accession}.log | reformat.sh in=stdin.fastq int=t verifypaired=t out1="${sample_ID}_R1.fastq.gz" out2="${sample_ID}_R2.fastq.gz" 2>>${accession}.log
done

# WARNING! This might not work if the R1 and R2 reads are different lengths... but I know that my data is okay.
```

Done! Next step is assembly and binning.

