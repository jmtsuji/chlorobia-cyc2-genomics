# ABOUT *cyc2* profile Hidden Markov Model development
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
conda create -y -n hmmer -c bioconda hmmer=3.2 clustalo=1.2.3 entrez-direct=11 seqtk=1.3
```
Use this environment via `conda activate hmmer`


## Get *cyc2* genes
```bash
work_dir="${github_repo_location}/Data_analysis_pipeline/04_HMM_development"
download_dir="${work_dir}/01_downloads"
guide_filepath="${download_dir}/cyc2_genome_info.tsv"
log_filepath="${download_dir}/cyc2_download.log"

mkdir -p ${download_dir}/final
cd ${download_dir}

# Get sequence info
# BE CAREFUL - some genomes have multiple cyc2's, separated by commas in the file. Must temporarily adjust the Internal Field Separator (IFS) - see below.
seq_names=($(cut -d $'\t' -f 1 ${guide_filepath} | tail -n +2))
seq_accessions_per_name=($(cut -d $'\t' -f 2 ${guide_filepath} | tail -n +2))

# Initialize output file
printf "" > ${download_dir}/cyc2_seqs_downloaded.faa
printf "" > ${log_filepath}

echo "[ $(date -u) ]: Downloading cyc2 from ${#seq_names[@]} genomes" | tee -a ${log_filepath}

for i in $(seq 1 ${#seq_names[@]}); do

    # Set counter to zero-ordered
    j=$((${i}-1))

    # Get variables
    seq_name=${seq_names[${j}]}

    # BUT for the seq accession, get all entries in case there are more than one (comma separated)
    OFS=${IFS}
    IFS=","
    seq_accessions=(${seq_accessions_per_name[${j}]})
    IFS=${OFS}

    echo "[ $(date -u) ]: Downloading ${#seq_accessions[@]} sequences for '${seq_name}'" | tee -a ${log_filepath}

    # Download
    for seq_accession in ${seq_accessions[@]}; do

        echo "[ $(date -u) ]: efetch -db sequences -id ${seq_accession} -format fasta" | tee -a ${log_filepath}
        efetch -db sequences -id ${seq_accession} -format fasta | seqtk seq -A | tail -n +2 > ${download_dir}/cyc2.tmp

        # Rename
        sed "1i >${seq_name}__${seq_accession}" ${download_dir}/cyc2.tmp >> ${download_dir}/cyc2_seqs_downloaded.faa
        rm ${download_dir}/cyc2.tmp

    done

done

echo "[ $(date -u) ]: Finished" | tee -a ${log_filepath}
```

Add the *Chlorobia* *cyc2* sequences from this study (already added to this folder)  
Note: in reality, these sequences were detected with an earlier version of the HMM, then the HMM was re-made to include these sequences later.
```bash
cd ${download_dir}
cat cyc2_seqs_metagenomes.faa cyc2_seqs_downloaded.faa > final/cyc2_seqs_all.faa
```

Make a separate file with only *Chlorobia*  
File `chlorobia_cyc2_seqs.list` contains the names of the *Chlorobia* sequences.
```bash
seqtk subseq final/cyc2_seqs_all.faa chlorobia_cyc2_seqs.list > final/cyc2_seqs_chlorobia.faa
```

Now ready to align and make the HMM.

## Align the genes
Using `clustalo` v1.2.3

```bash
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
```bash
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

