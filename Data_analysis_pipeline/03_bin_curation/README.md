# ABOUT bin dereplication and manual curation
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
conda create -y -n dRep -c bioconda -c r drep=2.0.5 r-dplyr
```
Use this environment via `conda activate dRep` (as shown below).



## Bin dereplication
Unique genome bins were determined using dRep, version 2.0.5.  
dRep was run for the lake metagenomes. The enrichment culture genomes were determined later to be distinct from any of the lake metagenomes with an additional dRep run, but then the lake metagenome dRep results were used, with the enrichment culture genomes appended to the output table (because of the timing of when data became available for the manuscript).

Prepare CheckM and bin information for dRep
```
#### Prepare genome bins and checkM summary info as input into dRep
# Set user variables
input_dir="${github_repo_location}/Data_analysis_pipeline/02_assembly_and_binning/lake_metagenomes"
output_dir="${github_repo_location}/Data_analysis_pipeline/03_bin_curation/01_dereplication"
threads=12

mkdir -p ${output_dir}/output ${output_dir}/input/bins

# TODO: assumes no samples have the same name
# Get all checkM output
echo "[ $(date -u) ]: finding CheckM output files"
completeness_files=($(find ${input_dir} -name completeness.tsv | grep "genomic_bins/checkm"))

# Get combined completeness/taxonomy header
echo "[ $(date -u) ]: getting CheckM table header"
sample_file=${completeness_files[0]}
join -a 1 -j 1 -t $'\t' ${sample_file} ${sample_file%/*}/taxonomy.tsv | head -n 1 > ${output_dir}/input/checkm_stats_all.tsv

# Join all completeness/taxonomy info
for c_file in ${completeness_files[@]}; do
    tax_file=${c_file%/*}/taxonomy.tsv
    echo "[ $(date -u) ]: Joining '${c_file}'"
    join -a 1 -j 1 -t $'\t' ${c_file} ${tax_file} | tail -n +2 >> ${output_dir}/input/checkm_stats_all.tsv
done

# Reduce and rename columns
echo "[ $(date -u) ]: simplifying CheckM table"
printf "genome\tcompleteness\tcontamination\tstrain_heterogeneity\ttaxonomy_contained\ttaxonomy_sister_lineage\n" > ${output_dir}/input/checkm_stats_reduced.tsv
cut -d $'\t' -f 1,12-14,18,19 ${checkm_filepath} | tail -n +2 | sed "s/\t/.fasta\t/" >> ${output_dir}/input/checkm_stats_reduced.tsv
# Also added .fasta to bin names to match filename, as required by dRep, using sed. Help from https://unix.stackexchange.com/a/36035, accessed 180406 at ~17:50 EDT.

# Convert to CSV and further reduce columns for dRep (just in case the program is picky)
cut -d $'\t' -f 1-3 ${output_dir}/input/checkm_stats_reduced.tsv | sed "s/\t/,/g" > ${output_dir}/input/checkm_stats_reduced.csv

# Make symlinks to genome bins
echo "[ $(date -u) ]: finding genome bins"
bins=($(cut -d $'\t' -f 1 ${output_dir}/input/checkm_stats_reduced.tsv | tail -n +2))

# TODO: make more robust. Is not terrible but could still have multiple search hits if someone named other files in this directory in the right way...
echo "[ $(date -u) ]: linking genome bins to dRep input folder"
for bin in ${bins[@]}; do
    sample_name=${bin%%.*}

    # Get bin path and filter by sample name
    bin_path=($(find ${input_dir} -type f -name ${bin} 2>/dev/null | grep "${sample_name}/genomic_bins/${bin}"))

    if [ ${#bin_path[@]} != 1 ]; then
        echo "ERROR: ${bin}: ${#bin_path[@]} matching bin paths (should be 1)."
        echo "Paths found (if any) are: ${bin_path[@]}"
        echo "Exiting..."
        exit 1
    fi

    echo "[ $(date -u) ]: ${bin}: ${bin_path[0]}"

    ln -s ${bin_path[0]} ${output_dir}/input/bins/${bin}
done
```

Dereplicate the bins
```
cd ${output_dir}/output
log_code=$(date '+%y%m%d_%H%M')
dRep dereplicate -p ${threads} -g ${output_dir}/input/bins/*.fasta --genomeInfo ${output_dir}/input/checkm_stats_reduced.csv ${output_dir}/output 2>&1 | tee dRep_${log_code}.log
```

Make a summary file of dRep output and CheckM input using R  
First, navigate to the right folder, then start R
```
cd ${output_dir}
R
```


```R
# Load libraries
library(dplyr)

# Read input tables
checkM <- read.table("input/checkm_stats_reduced.tsv", sep = "\t", header = T, stringsAsFactors = F)
Widb <- read.table("output/data_tables/Widb.csv", sep = ",", header = T, stringsAsFactors = F)
info <- read.table("output/data_tables/genomeInformation.csv", sep = ",", header = T, stringsAsFactors = F)

# Remove empty columns
Widb <- Widb[,-c(1,5,6,7,8,12,13,14,15)]

# Remove non-needed columns
info <- info[,-c(2,3)]

# Join tables
tbls_joined <- dplyr::left_join(Widb, info, by = "genome")
tbls_joined <- dplyr::left_join(tbls_joined, checkM, by = "genome")

# Reorder columns for clarity
tbls_joined <- tbls_joined[,c(2,3,1,4:ncol(tbls_joined))]

# Write output
write.table(tbls_joined, file = "output/dRep_summary.csv", sep = ",", row.names = F, col.names = T, quote = F)
quit(save = "no")
```
The output table here was cleaned up to produce `Table_S1.csv`, a copy of which is included in this folder.


## Manual cleanup of selected bins
All genome bins were imported into Anvi'o for subsequent manual cleanup of specific bins of interest. Bins of interest were seleted after the bin dereplication step based on taxonomy and CheckM stats (see manuscript).

### Lake metagenomes
Genome bins were imported into anvi'o using [atlas-to-anvi.sh](https://github.com/jmtsuji/atlas-extensions), version ____.

To install:
```
# TODO
```

To import the bins:
```
# Set user variables
input_dir="${github_repo_location}/Data_analysis_pipeline/03_bin_curation/01_dereplication"
output_dir="${github_repo_location}/Data_analysis_pipeline/03_bin_curation/02_anvio"
threads=12


```

To analyze a bin of interest (e.g., `___`):
```
# TODO
```

Once finished, export the genome bin sequences using
```
# TODO
```


## Contig ordering via mauve
Contigs for the manually curated genomes were ordered based on the reference genome sequence of *Chlorobium luteolum*. The re-ordering does not actually affect the downstream analyses of this paper but is good practice in some cases (e.g., for comparing gene arrangement).

Installed the mauve contig mover (a bit difficult):
```
# TODO
```

Got the curated genome bins (manually -- see instructions below)
```
# These should be in ${github_repo_location}/Data_analysis_pipeline/03_bin_curation/02_anvio/summarize if you made them youself.

# Place the genome bins into the folder ${github_repo_location}/Data_analysis_pipeline/03_bin_curation/03_curated_bins/unordered
# You'll have to create this folder
```

Downloaded *Chl. luteolum* genome
```
download_dir=""${github_repo_location}/Data_analysis_pipeline/03_bin_curation/03_curated_bins/reference"
luteolum_accession=

# TODO
```

Ordered the contigs
```
# TODO
```

Now, the genome bins are done!

## Taxonomic classification via GTDBTk
All genome bins were taxonomically classified using the Genome Tree Database Toolkit.

Downloaded the database
```
download_dir="${github_repo_location}/Data_analysis_pipeline/03_bin_curation/05_gtdbtk_taxonomy/database"

# TODO
```

Gathered together the relevant genome bins
```
destination_dir="${github_repo_location}/Data_analysis_pipeline/03_bin_curation/04_all_bin_summary"
# Need to change extension to ".fa" to work with one of the downstream scripts

# From the lake metagenomes

# From the enrichment culture metagenomes
```

Ran the classifier (needs ~100 GB RAM!)
```
# Set user variables
input_dir="${github_repo_location}/Data_analysis_pipeline/03_bin_curation/04_all_bin_summary"
output_dir="${github_repo_location}/Data_analysis_pipeline/03_bin_curation/05_gtdbtk_taxonomy"
threads=12

gtdbtk classify_wf --genome_dir ${input_dir} --out_dir ${output_dir} -x fa --min_perc_aa 0 --prefix ELA111314_dRep_gtdbtk --cpus ${threads}
```
The output files `ELA111314_dRep_gtdbtk.ar122.summary.tsv` and `ELA111314_dRep_gtdbtk.bac120.summary.tsv` will be used downstream in Figure 3.

Note: in the actual analysis, GTDBTk was run separately for the lake vs. enrichment culture metagenomes, and the output files were then concatenated. However, running them together is simpler and I think should not change the results much at all.

## Custom gene scanning via MetAnnotate
Several genes of interest were searched for in the entire binned genome set:
- *rpoB* - housekeeping gene (RNA polymerase gene) useful as a taxonomic marker
- *cyc2* - candidate gene marker for iron oxidation used in this study
- *dsrA* - a gene marker for sulfide oxidation or dissimilatory sulfate reduction

HMMs for *rpoB* and *dsrA* were downloaded from the FunGene website. The *cyc2* HMM was developed in this study and is available in Supplementary File S2 (`cyc2_all_references.hmm`).

Predicted amino acid sequences of all genome bins using prokka `v1.13.3`
```
input_dir="${github_repo_location}/Data_analysis_pipeline/03_bin_curation/04_all_bin_summary" # set up in GTDBTk step above
output_dir="${github_repo_location}/Data_analysis_pipeline/03_bin_curation/04_all_bin_summary/prokka"
threads=12

mkdir -p ${output_dir}
genome_bins=($(find ${input_dir} -maxdepth 1 -iname "*.fa" | sort -h))

for genome in ${genome_bins[@]}; do

    genome_base=${genome%.fa}
    genome_base=${genome_base##*/}

    # Clean up the name to remove any odd characters, to make more compatible downstream e.g., with MetAnnotate
    genome_base=$(echo ${genome_base} | awk '{ gsub("[^A-Za-z0-9]", "_"); print }')

    echo ${genome_base}
    prokka --outdir ${work_dir}/prokka/${genome_base} --prefix ${genome_base} --locustag ${genome_base}_ --cpus ${threads} ${genome} 2&>/dev/null
    # N.B., might have some problems with Archaea

done
```

Get HMMs (partially MANUAL)
```
hmm_dir="${github_repo_location}/Data_analysis_pipeline/03_bin_curation/06_metannotate/hmm_files"
threads=12

# 1. Download the rpoB and dsrA HMMs from FunGene and place them in the above directory!

# 2. Then make the cyc2 HMM or use the one in Supplementary File S2 (cyc2_all_references.hmm)

# To make: start with the amino acid alignment in Supplementary File S3 (cyc2_Fig1_unmasked.faa). Move to the folder specified below.
input_file="${github_repo_location}/Data_analysis_pipeline/03_bin_curation/06_metannotate/hmm_files_rough_work/cyc2_Fig1_unmasked.faa"
output_dir="${github_repo_location}/Data_analysis_pipeline/03_bin_curation/06_metannotate/hmm_files_rough_work"

file_base=${input_file%.*}

# De-align sequences (just to show what I did)
seqtk seq -A ${file} | awk '{ if ($0 !~ /^>/) { gsub("-", ""); } print }' > ${file_base}_unaligned.faa

# Re-align sequences
cd ${output_dir}
clustalo --full --percent-id --wrap=60 --distmat-out=${file_base}_aligned_distmat.txt -i ${file_base}_unaligned.faa -o ${file_base}_aligned.faa --threads ${threads} -v > ${file_base}_aligned.log

# Build HMM
echo "hmmbuild -O ${file_base}.stk -n ${file_base} --cpu ${threads} ${file_base}.hmm ${file_base}_unaligned.faa > ${file_base}_hmmbuild.log" | tee ${file_base}_hmmbuild.log
hmmbuild -O ${file_base}.stk -n ${file_base} --cpu ${threads} ${file_base}.hmm ${input_file} >> ${file_base}_hmmbuild.log

done

# TODO - add sequence collection information (e.g., how to download)
```

Downloaded the RefSeq protein database  
**SKIP** this step if you've already downloaded the database before
```
# User variables
orf_dir="${github_repo_location}/Data_analysis_pipeline/03_bin_curation/04_all_bin_summary/prokka"
hmm_dir="${github_repo_location}/Data_analysis_pipeline/03_bin_curation/06_metannotate/hmm_files"
output_dir="${github_repo_location}/Data_analysis_pipeline/03_bin_curation/06_metannotate/output_contigs"
download_dir="${github_repo_location}/Data_analysis_pipeline/03_bin_curation/06_metannotate/database"

mkdir -p ${download_dir} ${output_dir}

# Enter the docker container (will auto-download if you don't have it installed)
enter-metannotate ${download_dir} ${orf_dir} ${hmm_dir} ${output_dir}

# Download the database (can take a while)
cd $METANNOTATE_DIR
sudo chown linuxbrew ../databases
bash refseq_installation.sh /home/linuxbrew/databases
sudo chown -R root:root ../databases
# The 'chown' commands temporarily make the output folder belong to the linuxbrew user inside the Docker container so that the user can run the Docker commands. Files are given back to you at the end.

# Leave the container
exit
```

Ran MetAnnotate on the predicted amino acid sequences
```
# User variables
# NOTE: CHANGE the database_dir path if you've already downloaded Refseq.fa somewhere else on your machine and want to use that instead
orf_dir="${github_repo_location}/Data_analysis_pipeline/03_bin_curation/04_all_bin_summary/prokka"
hmm_dir="${github_repo_location}/Data_analysis_pipeline/03_bin_curation/06_metannotate/hmm_files"
output_dir="${github_repo_location}/Data_analysis_pipeline/03_bin_curation/06_metannotate/output_contigs"
download_dir="${github_repo_location}/Data_analysis_pipeline/03_bin_curation/06_metannotate/database"

# Enter the docker container (will auto-download if you don't have it installed)
enter-metannotate ${download_dir} ${orf_dir} ${hmm_dir} ${output_dir}

threads=12

ref_UID=$(stat -c "%u" /home/linuxbrew/output)
sudo chown -R linuxbrew /home/linuxbrew/output
echo ${threads} > MetAnnotate/concurrency.txt # set the number of threads
metannotate-wrapper-docker sequence orf_files hmm_files 2>&1 | tee output/metannotate_wrapper_docker.log
sudo chown -R $ref_UID /home/linuxbrew/output
# The 'chown' commands temporarily make the output folder belong to the linuxbrew user inside the Docker container so that the user can run the Docker commands. Files are given back to you at the end.

# Leave the docker container
exit
```

The output file `all_annotations`[random_code]`.tsv` shows the hits of the HMMs on all genomes, along with the e-value score. This is used in the code for Figure 3, and a copy of the file is the Figure 3 folder.

## Relative abundance profiling using read mapping
### First had to install (only for first use)
```
## Get the git repo
work_dir="${github_repo_location}/Data_analysis_pipeline/03_bin_curation/07_read_mapping"
mkdir -p ${work_dir} && cd ${work_dir}

git clone https://github.com/jmtsuji/atlas-extensions.git
cd ${work_dir}/atlas-extensions
git checkout perfectmode

# TODO - checkout specific version # instead.

## Make the conda env and print the versions
conda create -y -n bin_mapping_stats -c conda-forge -c bioconda -c r samtools bbmap r r-plyr r-dplyr
conda install -n bin_mapping_stats -c r r-getopt # neglected to do this the first time around
conda list -n bin_mapping_stats > ../conda_package_versions_bin_mapping_stats.txt

## Finish install
conda activate bin_mapping_stats
env_path=$(which samtools)
env_path=${env_path%/*}
cp calculate_bin_abundance_in_metagenome.sh calculate_coverage_stats.R aggregate_mapping_stats.R ${env_path}
```

Then run the profiler
```
# TODO - not yet cleaned up

cd ${work_dir}
calculate_bin_abundance_in_metagenome.sh input/dereplicated_genomes input/metagenome_QC_reads output ${threads} ${memory} 2>&1 | tee calculate_bin_abundance_in_metagenome.log
# Started 180927 at ~2:52 AM EDT

# Ran into issue where the script failed if there were 0 mapped reads. Corrected for this in commit 5b14f12
# Used this commit to continue. Did not want to repeat metagenome read counting (code not changed here), so skipped this step
# by creating a copy of the script called 'calculate_bin_abundance_in_metagenome_part2.sh' with some lines commented out.
# Delete mapping and coverage folders, and the coverage and mapping summary files.
# Then started again:
cd ${out_dir}
./calculate_bin_abundance_in_metagenome_part2.sh input/dereplicated_genomes input/metagenome_QC_reads output ${threads} ${memory} 2>&1 | tee calculate_bin_abundance_in_metagenome.log.2
# Started 180927, ~3:20 AM EDT
```

## Relative abundance profiling of raw reads via MetAnnotate
This was done as a cross-comparison to the above code.

The same three HMMs as mentioned earlier for MetAnnotate were used, along with one more for reference: *bchL* - pigment synthesis gene. This HMM was also downloaded from Fun
```
# TODO
```

### Done!
This is the end of the main heavy-lifting data processing work for this paper. Figures were generated based off this dataset using code found in each figure folder.

