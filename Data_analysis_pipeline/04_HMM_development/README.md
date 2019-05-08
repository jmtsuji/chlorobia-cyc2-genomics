# ABOUT *cyc2* profile Hidden Markov Model development
Copyright Jackson M. Tsuji, Neufeld Research Group, 2019
Part of the larger *IISD-ELA Chlorobia cyc2 project*.

All code here is to be run in a Bash terminal of a unix system (e.g., Ubuntu 16 or 18) unless indicated otherwise.

## Define where you downloaded the Github repo:
```
github_repo_location="/Analysis/jmtsuji/Chlorobia_cyc2_code"
```

## Software prerequisites
- miniconda (miniconda3 preferred)


## Conda environment with all needed dependencies:
```
conda create -y -n hmmer -c bioconda hmmer=3.2 clustalo=1.2.3 entrez-direct=11 seqtk=1.3-r106

# Also manually install a helpful renaming script
cd /tmp
wget -O - https://github.com/jmtsuji/basic-sequence-analysis/archive/v1.2.0.tar.gz | tar -xzf -
# Creates directory 'basic-sequence-analysis-1.2.0'
cd basic-sequence-analysis-1.2.0
cp basic-sequence-analysis-version text_find_and_replace.sh ${CONDA_PREFIX}/bin
cd ..
rm -rf basic-sequence-analysis-1.2.0
```
Use this environment via `conda activate hmmer`


## Download genomes and identify *cyc2* genes
```
work_dir="${github_repo_location}/Data_analysis_pipeline/04_HMM_development"
download_dir="${work_dir}/01_downloads"
guide_filepath="${download_dir}/cyc2_genome_info.tsv"
log_filepath="${download_dir}/cyc2_download.log"

mkdir -p ${download_dir}/final
cd ${download_dir}

# Get sequence info
seq_names=($(cut -d $'\t' -f 1 ${guide_filepath} | tail -n +2))
seq_accessions=($(cut -d $'\t' -f 2 ${guide_filepath} | tail -n +2))

# Initialize output file
printf "" > ${download_dir}/cyc2_seqs_downloaded_raw.faa
printf "" > ${log_filepath}

echo "[ $(date -u) ]: Downloading ${#seq_names[@]} sequences" | tee -a ${log_filepath}

for i in $(seq 1 ${#seq_names[@]}); do
# Set counter to zero-ordered
j=$((${i}-1))

# Get variables
seq_name=${seq_names[${j}]}
seq_accession=${seq_accessions[${j}]}

# Download
echo "[ $(date -u) ]: Downloading '${seq_accession}' for '${seq_name}'" | tee -a ${log_filepath}
echo "[ $(date -u) ]: efetch -db sequences -id ${seq_accession} -format fasta >> ${download_dir}/cyc2_seqs_downloaded_raw.faa" | tee -a ${log_filepath}
efetch -db sequences -id ${seq_accession} -format fasta >> ${download_dir}/cyc2_seqs_downloaded_raw.faa

done

echo "[ $(date -u) ]: Finished" | tee -a ${log_filepath}
```

Rename the cyc2 seqences for clarity
```
# Get rid of comments
cd ${download_dir}
seqtk seq -C cyc2_seqs_downloaded_raw.faa > cyc2_seqs_downloaded_no_comments.faa.tmp

# Make renaming guide file
printf "old_names\n" > renaming_guide.tsv.tmp.1
cut -d $'\t' -f 2 ${guide_filepath} | tail -n +2 >> renaming_guide.tsv.tmp.1

printf "new_names\n" > renaming_guide.tsv.tmp.2
cut -d $'\t' -f 1 ${guide_filepath} | tail -n +2 >> renaming_guide.tsv.tmp.2

paste renaming_guide.tsv.tmp.1 renaming_guide.tsv.tmp.2 > renaming_guide.tsv

# Find and replace sequence headers
text_find_and_replace.sh renaming_guide.tsv cyc2_seqs_downloaded_no_comments.faa.tmp cyc2_seqs_downloaded_renamed.faa 2>&1 | tee text_find_and_replace.log
rm cyc2_seqs_downloaded_no_comments.faa.tmp renaming_guide.tsv.tmp.1 renaming_guide.tsv.tmp.2
```

Add the *Chlorobia* *cyc2* sequences from this study (already added to this folder)  
Note: in reality, these sequences were detected with an earlier version of the HMM, then the HMM was re-made to include these sequences later.
```
cd ${download_dir}
cat cyc2_seqs_metagenomes.faa cyc2_seqs_downloaded_renamed.faa > final/cyc2_seqs_all.faa
```

Make a separate file with only *Chlorobia* (manually)  
File `chlorobia_cyc2_seqs.list` contains the names of the *Chlorobia* sequences.
```
seqtk subseq final/cyc2_seqs_all.faa chlorobia_cyc2_seqs.list > final/cyc2_seqs_chlorobia.faa
```

Now ready to align and make the HMM.

## Align the genes
Using `clustalo` v1.2.3

```
work_dir="${github_repo_location}/Data_analysis_pipeline/04_HMM_development"
input_dir="${work_dir}/01_downloads/final"
input_filepaths=(${input_dir}/cyc2_seqs_all.faa ${input_dir}/cyc2_seqs_chlorobia.faa) # HARD-CODED
alignment_dir="${work_dir}/02_alignment"
threads=4

mkdir -p ${alignment_dir}
cd ${alignment_dir}

for input_file in ${input_filepaths[@]}; do

base_name=${input_file%.faa}
base_name=${base_name##*/}
echo "[ $(date -u) ]: Aligning sequence '${base_name}.faa'"

clustalo --in ${input_file} --out ${base_name}_aligned.faa --distmat-out=${base_name}_distmat.txt --full --percent-id --threads=${threads} --verbose 2>&1 | tee ${base_name}_alignment.log

done
```

## Create the HMMs
```
work_dir="${github_repo_location}/Data_analysis_pipeline/04_HMM_development"
input_dir="${work_dir}/02_alignment"
input_filepaths=(${input_dir}/cyc2_seqs_all_aligned.faa ${input_dir}/cyc2_seqs_chlorobia_aligned.faa) # HARD-CODED
hmm_names=(cyc2_all_references cyc2_Chlorobia) # HARD-CODED to mathc with above
hmm_dir="${work_dir}/03_HMMs"

mkdir -p ${hmm_dir}
cd ${hmm_dir}

for i in $(seq 1 ${#input_filepaths[@]}); do
# Set counter to zero-ordered
j=$((${i}-1))

# Get variables
input_filepath=${input_filepaths[${j}]}
hmm_name=${hmm_names[${j}]}

# Build HMM
echo "[ $(date -u) ]: Building HMM '${hmm_name}' for alignment '${input_filepath##*/}'"
hmmbuild -n ${hmm_name} -o ${hmm_name}.log -O ${input_filepath}.mod ${hmm_name}.hmm ${input_filepath}
done

```

HMM finished!

