# ABOUT calculating bin environmental abundance from dRep bins
ELA 2011-2014 metagenomes
Jackson M. Tsuji, 2018
Background: wrote code in atlas-extensions repo to be able to calculate environmental abundances in an automated way. Trying this out for the ELA dataset.

# 1. Ran for the dRep'ed unrefined genome bins in the ELA 2011-2014 dataset, plus refined bins and enrichments
Working on mellea
Using commit 5b14f12 in jmtsuji/atlas-extensions (branch 'perfectmode')
BUT started on commit f600173 (for metagenome read counting; exact same code as the later commit - see notes below)

```
work_dir="/Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full"
out_dir=${work_dir}/post-analysis/11_bin_env_abundance/vs2_dRep_and_enrichments
threads=20
memory=100

dRep_bins_tarball=${work_dir}/post-analysis/05_dRep/03b_all_assemblies_updated_L442/output/dereplicated_genomes.tar.bz2

# Grab dRep'ed genomes and change extension to .fa
mkdir -p ${out_dir}/input
bin_out_dir=${out_dir}/input/dereplicated_genomes
tar -xvjf ${dRep_bins_tarball} -C ${out_dir}/input

# Renamed to "*.fa"
genomes=($(find ${out_dir}/input/dereplicated_genomes -iname "*.fasta"))
for genome in ${genomes[@]}; do
echo ${genome##*/}
mv ${genome} ${genome%.*}.fa
done

# Moved unrefined bins to a separate folder
mkdir -p ${out_dir}/bins_old
mv ${bin_out_dir}/CA-L227-2013.22.fa ${bin_out_dir}/CA-L227-2013.56.fa \
	${bin_out_dir}/CA-L227-2013.55.fa ${bin_out_dir}/CA-L442.64.fa ${out_dir}/bins_old

# Grabbed refined genomes
# From ELA
ELA_bin_dir="/Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/post-analysis/06_refined_bins/01_Chlorobi/refined_FINAL"
ELA_bins=($(find ${ELA_bin_dir} -iname "*.fa"))
for bin in ${ELA_bins[@]}; do
echo ${bin##*/}
ln ${bin} ${bin_out_dir}/${bin##*/}
done
# From enrichments
enr_bin_dir="/Analysis/jmtsuji/Documents/02_Active_workfolders/40_180721_GSB_genome_seq_analysis/post-analysis/04_refined_bins/01_Chlorobi"
enr_bins=($(find ${enr_bin_dir} -iname "*.fa"))
for bin in ${enr_bins[@]}; do
echo ${bin##*/}
ln ${bin} ${bin_out_dir}/${bin##*/}
done

# Grab metagenomes (hard links)
mkdir -p ${out_dir}/input/metagenome_QC_reads
metagenomes=($(find ${work_dir} -iname "*QC*fastq.gz" | grep -v "post-analysis" | grep -v "coassembly"))
for metagenome in ${metagenomes[@]}; do
	ln ${metagenome} ${out_dir}/input/metagenome_QC_reads/${metagenome##*/}
	echo ${metagenome##*/}
done
```

Install code
```
#### INSTALL - only for the first time
## Get the git repo
# git clone https://github.com/jmtsuji/atlas-extensions.git
# cd ${out_dir}/atlas-extensions
# git checkout perfectmode

## Only for first time: make the conda env and print the versions
# conda create -y -n bin_mapping_stats -c conda-forge -c bioconda -c r samtools bbmap r r-plyr r-dplyr
# conda install -n bin_mapping_stats -c r r-getopt # neglected to do this the first time around
# conda list -n bin_mapping_stats > ../conda_package_versions_bin_mapping_stats.txt

## As needed
# source activate bin_mapping_stats

## Finish install
# env_path=$(which samtools)
# env_path=${env_path%/*}
# cp calculate_bin_abundance_in_metagenome.sh calculate_coverage_stats.R aggregate_mapping_stats.R ${env_path}
```

Run the environmental mapping
```
# Run!
cd ${out_dir}
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


# 2. Ran GTDB classifier on the genome bin set
Used classifier version 'GTDB-Tk v0.1.3' with database release 86 (Sept. 2018)
Setup script for conda env is in a Gitub repo of mine: https://github.com/jmtsuji/gtdbtk_testing (I'm almost sure I used this method)

Command (from log file):
```
gtdbtk classify_wf --genome_dir /Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/post-analysis/11_bin_env_abundance/vs2_dRep_and_enrichments/input/dereplicated_genomes --out_dir /Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/post-analysis/11_bin_env_abundance/vs2_dRep_and_enrichments/input/dereplicated_genomes/gtdbtk -x fa --min_perc_aa 0 --prefix ELA111314_dRep_gtdbtk --cpus 40
```

Output at `/Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/post-analysis/11_bin_env_abundance/vs2_dRep_and_enrichments/input/gtdbtk`

# 3. Re-made the cyc2 HMM with several variants for additional analysis

# Side note: making new HMMs
folder: `cyc2_hmm` on my laptop, inside the env abundance analysis folder: `/home/jmtsuji/Research_General/PhD/04b_Metagenome_resequencing_F2015/10_ATLAS_re_analysis/09_env_abundance_of_bins/vs2_dRep_scripted/input/cyc2_hmm`

Software versions:
- Clustalo 1.2.3
- hmmbuild 3.1b2
- seqtk 1.3-r106

Inputs (from the related cyc2 comparison analysis):
- 3GSB -- using just C ferro, C phaeo, C luteolum. Mimicking the one I made a year ago (should hopefully be identical)
- 7GSB -- including the four cyc2 sequences from the most recent study
- ren1_noLepto_3GSB -- Includes all reference Fe oxidizers (excluding the distant Leptospirillum) but OMITS the 4 from this study
- ren1_noLepto -- already aligned. Includes all reference Fe oxidizers (excluding the distant Leptospirillum) plus the 4 from this study


```
files=(cyc2_comparison_collection_3GSB.faa cyc2_comparison_collection_7GSB.faa cyc2_whole_sequence_set_ren1_noLepto_3GSB.faa)

for file in ${files[@]}; do

file_base=${file%.*}
echo ${file_base}

# De-align sequence
seqtk seq -A ${file} | awk '{ if ($0 !~ /^>/) { gsub("-", ""); } print }' > ${file_base}_unaligned.faa

# Re-align sequence
clustalo --full --percent-id --wrap=60 --distmat-out=${file_base}_aligned_distmat.txt -i ${file_base}_unaligned.faa -o ${file_base}_aligned.faa --threads 4 -v > ${file_base}_aligned.log

done

files=(cyc2_comparison_collection_3GSB_aligned.faa cyc2_comparison_collection_7GSB_aligned.faa cyc2_whole_sequence_set_ren1_noLepto_3GSB_aligned.faa cyc2_whole_sequence_set_ren1_noLepto_aligned.faa)
names=(cyc2_3GSB cyc2_7GSB cyc2_reference_genomes cyc2_reference_genomes_plus4GSB)

for i in $(seq 1 ${#files[@]}); do
j=$((${i}-1))
file=${files[${j}]}
name=${names[${j}]}
echo "${name}: ${file}"

# Build HMM
echo "hmmbuild -O ${name}.stk -n ${name} --cpu 4 ${name}.hmm ${file} > ${name}_hmmbuild.log" | tee ${name}_hmmbuild.log
hmmbuild -O ${name}.stk -n ${name} --cpu 4 ${name}.hmm ${file} >> ${name}_hmmbuild.log

done

conda list > conda_env.log

```

Copied these HMMs to the server at: `/Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/post-analysis/11_bin_env_abundance/vs2_dRep_and_enrichments/input/cyc2_hmm`


# 4. Added two extra GSB genome bins to the overall collection
These are dirtier bins from the L227 assemblies that did not pass the dRep threshold, but I think they are important to include in the environmental analysis.

Made folder `/Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/post-analysis/11_bin_env_abundance/vs2_dRep_and_enrichments/input/dereplicated_genomes/extra_Chlorobi_genomes_dirty`
Added `CA-L227-2014.106.fasta  CA-L227-2014.92.fasta` from the L227 coassembly

A couple notes:
- 106: this one seems very similar to 2013-71. Both are ~90% complete but ~70% contaminated. Have two dsrA, I think.
- 92: this one does not appear in 2013. Has two cyc2's (maybe the second is unreliable?). 77% complete, 27% contaminated.

## More data on these bins
```

bin_id	checkm_bin_taxonomy_contained	checkm_bin_taxonomy_sister_lineage	checkm_bin_number_unique_markers	checkm_bin_completeness	checkm_bin_contamination	L227.2013.6m_count	L227.2013.8m_count	L227.2014.6m_count	L227.2014.8m_count
CA-L227-2014.92	k__Bacteria;p__Chlorobi;c__Chlorobia;o__Chlorobiales;f__Chlorobiaceae;g__Pelodictyon	s__Pelodictyon_phaeoclathratiforme	10	77.45	27.64	5.5	26.5	124	15
CA-L227-2014.92	k__Bacteria;p__Chlorobi;c__Chlorobia;o__Chlorobiales;f__Chlorobiaceae;g__Pelodictyon	s__Pelodictyon_phaeoclathratiforme	10	77.45	27.64	2.5	15.5	135.5	20

bin_id	checkm_bin_taxonomy_contained	checkm_bin_taxonomy_sister_lineage	checkm_bin_number_unique_markers	checkm_bin_completeness	checkm_bin_contamination	L227.2013.6m_count	L227.2013.8m_count	L227.2014.6m_count	L227.2014.8m_count
CA-L227-2013.71	k__Bacteria;p__Chlorobi;c__Chlorobia;o__Chlorobiales;f__Chlorobiaceae;g__Pelodictyon	s__Pelodictyon_phaeoclathratiforme	15	91.38	69.44	56.5	127.5	404.5	62
CA-L227-2014.106	k__Bacteria;p__Chlorobi;c__Chlorobia;o__Chlorobiales;f__Chlorobiaceae;g__Pelodictyon	s__Pelodictyon_phaeoclathratiforme	13	91.14	67.16	51.5	127	356	51

```

Other Chlorobi bins besides the four main ones will also be included in the analysis -- these got through dRep but were not good enough for what I was looking for with the rest of the paper. E.g., CA-L442 bin 74.


# 5. Realized I should expand the read mapping approach to include more genomes and more metagenomes
- More genomes: the two extras added above in step 4
- More metagenomes: the enrichment culture metagenomes

## Genomes already added above

## Added the enrichment metagenomes
Put into subdirectory at `/Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/post-analysis/11_bin_env_abundance/vs2_dRep_and_enrichments/input/metagenome_QC_reads/enrichments`

Made via hard links:
```
dest_dir="/Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/post-analysis/11_bin_env_abundance/vs2_dRep_and_enrichments/input/metagenome_QC_reads/enrichments"

source_dir="/Analysis/jmtsuji/Documents/02_Active_workfolders/40_180721_GSB_genome_seq_analysis/L227_S_6D/sequence_quality_control"
files=(L227_S_6D_QC_R1.fastq.gz L227_S_6D_QC_R2.fastq.gz L227_S_6D_QC_se.fastq.gz)

for file in ${files[@]}; do
ln ${source_dir}/${file} ${dest_dir}/${file}
done


source_dir="/Analysis/jmtsuji/Documents/02_Active_workfolders/40_180721_GSB_genome_seq_analysis/L304_S_6D/sequence_quality_control"
files=(L304_S_6D_QC_R1.fastq.gz L304_S_6D_QC_R2.fastq.gz L304_S_6D_QC_se.fastq.gz)

for file in ${files[@]}; do
ln ${source_dir}/${file} ${dest_dir}/${file}
done

```

Each analysis runs independently, so in theory I could just run these additional analysis separated and combine into the same table as the main run

## Ran the analysis for the two new genomes against the ELA metagenomes
```
source activate bin_mapping_stats
cd /Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/post-analysis/11_bin_env_abundance/vs2_dRep_and_enrichments

threads=10
memory=80

# Looks like still on commit 5b14f12, like the run above (more or less)
calculate_bin_abundance_in_metagenome.sh input/dereplicated_genomes/extra_Chlorobi_genomes_dirty input/metagenome_QC_reads output/calculate_bin_abundance_in_metagenome_part2 ${threads} ${memory} 2>&1 | tee calculate_bin_abundance_in_metagenome_part2_extraChlorobi_181206.log
# For this to work, I had to temporarily move the enrichment metagenomes from `input/metagenome_QC_reads/enrichments` to `input/metagenome_QC_reads`. Bug in script.
```

## Added enrichment bins to the dRep pool:
# Re-ran dRep to include the two new metagenomes of enrichment cultures
```
source_dir="/Analysis/jmtsuji/Documents/02_Active_workfolders/40_180721_GSB_genome_seq_analysis"
output_dir="/Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/post-analysis/05_dRep/04_with_J17_2GSB_enrichments"
threads=12

# Make a copy of the old input dir for reference
cp -r /Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/post-analysis/05_dRep/03b_all_assemblies_updated_L442/input ${output_dir}

# Fix the broken symlinks
cd ${output_dir}/input/bins
names=($(find $PWD -iname "*.fasta" | sort -h))
links_old=($(find $PWD -iname "*.fasta" | sort -h | xargs readlink))

## Example old link
# /Hippodrome/jmtsuji/180123_ELA111314_atlas_r1.0.22_plus_full/L227-2013-6m/genomic_bins/L227-2013-6m.025.fasta
## What it should be
# /Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/L227-2013-6m/genomic_bins/L227-2013-6m.025.fasta

printf "" > link_repair.tsv
for i in $(seq 1 ${#names[@]}); do
j=$((${i}-1))
destination=${names[${j}]}
link_old=${links_old[${j}]}
link_new=${link_old##/Hippodrome/jmtsuji/}
link_new="/Analysis/jmtsuji/Hippodrome_legacy/${link_new}"
printf "${destination}\t${link_new}\n" >> link_repair.tsv
done

links_new=($(cut -d $'\t' -f 2 link_repair.tsv))
for i in $(seq 1 ${#names[@]}); do
j=$((${i}-1))
destination=${names[${j}]}
link_new=${links_new[${j}]}

rm ${destination}
ln -s ${link_new} ${destination}

done
# wow! it looks like it actually worked.

# Get checkm info for the new new runs
# Combine checkM output from the latest runs
completeness_files=(${source_dir}/L304_S_6D/genomic_bins/checkm/completeness.tsv ${source_dir}/L227_S_6D/genomic_bins/checkm/completeness.tsv)
head -n 1 ${output_dir}/input/checkm_stats_all.tsv > ${output_dir}/input/checkm_stats_all_enrichments.tsv
for c_file in ${completeness_files[@]}; do
tax_file=${c_file%/*}/taxonomy.tsv
echo $c_file
join -a 1 -j 1 -t $'\t' ${c_file} ${tax_file} | tail -n +2 >> ${output_dir}/input/checkm_stats_all_enrichments.tsv
done

checkm_filepath=${output_dir}/input/checkm_stats_all_enrichments.tsv

# Reduce and rename columns
printf "genome\tcompleteness\tcontamination\tstrain_heterogeneity\ttaxonomy_contained\ttaxonomy_sister_lineage\n" > ${output_dir}/input/checkm_stats_reduced_enrichments.tsv
cut -d $'\t' -f 1,12-14,18,19 ${checkm_filepath} | tail -n +2 | sed "s/\t/.fasta\t/" >> ${output_dir}/input/checkm_stats_reduced_enrichments.tsv
# Also added .fasta to bin names to match filename, as required by dRep, using sed. Help from https://unix.stackexchange.com/a/36035, accessed 180406 at ~17:50 EDT.

# Convert to CSV and further reduce columns for dRep (just in case the program is picky)
cut -d $'\t' -f 1-3 ${output_dir}/input/checkm_stats_reduced_enrichments.tsv | sed "s/\t/,/g" > ${output_dir}/input/checkm_stats_reduced_enrichments.csv

# Added checkm info to the main table
cp checkm_stats_reduced.csv checkm_stats_reduced_ORIGINAL_BACKUP.csv # made a backup
tail -n +2 ${output_dir}/input/checkm_stats_reduced_enrichments.csv >> checkm_stats_reduced.csv

# Get symlinks to genome bins -- just did this manually
ln -s ${source_dir}/L304_S_6D/genomic_bins/L304_S_6D.001.fasta ${output_dir}/input/bins/L304_S_6D.001.fasta
ln -s ${source_dir}/L227_S_6D/genomic_bins/L227_S_6D.001.fasta ${output_dir}/input/bins/L227_S_6D.001.fasta
ln -s ${source_dir}/L227_S_6D/genomic_bins/L227_S_6D.002.fasta ${output_dir}/input/bins/L227_S_6D.002.fasta
ln -s ${source_dir}/L227_S_6D/genomic_bins/L227_S_6D.003.fasta ${output_dir}/input/bins/L227_S_6D.003.fasta

#### Run dRep
mkdir -p ${output_dir}/output
cd ${output_dir}/output
source activate dRep
conda list > conda_env.log
# Using dRep v2.0.5
log_code=$(date '+%y%m%d_%H%M')
dRep dereplicate -p ${threads} -g ${output_dir}/input/bins/*.fasta --genomeInfo ${output_dir}/input/checkm_stats_reduced.csv ${output_dir}/output 2>&1 | tee dRep_${log_code}.log

### Make a summary file of dRep output and CheckM input using R
cd ${output_dir}

R
library(dplyr)
checkM <- read.table("input/checkm_stats_reduced.tsv", sep = "\t", header = T, stringsAsFactors = F)
Widb <- read.table("output/data_tables/Widb.csv", sep = ",", header = T, stringsAsFactors = F)
info <- read.table("output/data_tables/genomeInformation.csv", sep = ",", header = T, stringsAsFactors = F)

# Remove empty columns
Widb <- Widb[,-c(1,5,6,7,8,12,13,14,15)]

# Remove non-needed columns
info <- info[,-c(2,3)]

# Join
tbls_joined <- dplyr::left_join(Widb, info, by = "genome")
tbls_joined <- dplyr::left_join(tbls_joined, checkM, by = "genome")

# Reorder columns for clarity
tbls_joined <- tbls_joined[,c(2,3,1,4:ncol(tbls_joined))]

write.table(tbls_joined, file = "output/dRep_summary.csv", sep = ",", row.names = F, col.names = T, quote = F)
quit(save = "no")

```
Results were a bit baffling. The two new genomes were added in the final rep set, BUT the final rep set was ~50 genomes larger than the last time I ran the tool,
also on version 2.0.5. So odd. I need to look into this sometime to see what happened.
But for now, I will assume that the two new genomes (from the L227 enrichment) are unique compared to the other dRep'ed genomes and add them to the analysis.

## Ran the analysis for the two new enrichment genomes against the ELA metagenomes, like what I just did for the new GSB genomes above.
```
source activate bin_mapping_stats
cd /Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/post-analysis/11_bin_env_abundance/vs2_dRep_and_enrichments

threads=10
memory=80

# Looks like still on commit 5b14f12, like the run above (more or less)
# This time saved time by using the read counts already performed a few minutes ago above and just mapping reads, using the quick script already in this folder
./calculate_bin_abundance_in_metagenome_part2.sh input/dereplicated_genomes/extra_genomes_enrichments input/metagenome_QC_reads output/calculate_bin_abundance_in_metagenome_part3 ${threads} ${memory} 2>&1 | tee calculate_bin_abundance_in_metagenome_part3_enr_bins_181206_0413.log
# For this to work, I had to temporarily move the enrichment metagenomes from `input/metagenome_QC_reads/enrichments` to `input/metagenome_QC_reads`. Bug in script.

# Right after this started running (outside the screen), copied the metagenome read counts table
cp output/calculate_bin_abundance_in_metagenome_part2/details_stats/metagenome_read_counts.tsv output/calculate_bin_abundance_in_metagenome_part3/details_stats
statswrapper.sh input/dereplicated_genomes/extra_genomes_enrichments/*.fa format=5 > output/calculate_bin_abundance_in_metagenome_part3/detailed_stats/genome_stats.tsv
# This is not ideal but works...
```

## Then got read mapping information for all other genomes (excluding these four) against the enrichment metagenomes
Moved the four genomes temporarily out of the `input/dereplicated_genomes` folder to accomplish this.
```
source activate bin_mapping_stats
cd /Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/post-analysis/11_bin_env_abundance/vs2_dRep_and_enrichments

threads=20
memory=80

# Looks like still on commit 5b14f12, like the run above (more or less)
calculate_bin_abundance_in_metagenome.sh input/dereplicated_genomes input/metagenome_QC_reads/enrichments output/calculate_bin_abundance_in_metagenome_part4 ${threads} ${memory} 2>&1 | tee calculate_bin_abundance_in_metagenome_part4_enr_metagenomes_181206_0421.log

```


## Also classified the four new genomes via gtdbtk:
Temporarily moved into the same dir for ease of use
```
source activate gtdbtk

genome_dir="/Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/post-analysis/11_bin_env_abundance/vs2_dRep_and_enrichments/input/temp/new"
out_dir="/Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/post-analysis/11_bin_env_abundance/vs2_dRep_and_enrichments/output/gtdbtk_4new_181206"
mkdir -p ${out_dir}
conda list > ${out_dir}/conda_env.log

gtdbtk classify_wf --genome_dir ${genome_dir} --out_dir ${out_dir} -x fa --min_perc_aa 0 --prefix ELA_J17_enr_gtdbtk --cpus 20

```


### ALL DONE

Will combine on my laptop -- see later step


# 6. Identify genes in genome bins using MetAnnotate
Basic process:
- predict amino acids using prokka (to allow me to double check metannotate annotations if desired)
- run MetAnnotate

Software:
- prokka: `1.13.3`

## Part 1: predicted proteins
```
# Predicted proteins
work_dir="/Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/post-analysis/11_bin_env_abundance/vs2_dRep_and_enrichments/input/dereplicated_genomes"

mkdir -p ${work_dir}/prokka

# All genomes have extensions of .fa
genomes=($(find ${work_dir} -name "*fa" -type f | sort -h))

for genome in ${genomes[@]}; do

genome_base=${genome%.fa}
genome_base=${genome_base##*/}

# Clean up the name to remove any odd characters, to make more compatible downstream e.g., with MetAnnotate
genome_base=$(echo ${genome_base} | awk '{ gsub("[^A-Za-z0-9]", "_"); print }')

echo ${genome_base}
prokka --outdir ${work_dir}/prokka/${genome_base} --prefix ${genome_base} --locustag ${genome_base}_ --cpus 40 ${genome} 2&>/dev/null
# N.B., might have some problems with Archaea

done
```

### Added the two new genomes from the enrichments
```
work_dir="/Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/post-analysis/11_bin_env_abundance/vs2_dRep_and_enrichments/input/dereplicated_genomes"
genomes=($(find ${work_dir}/extra_genomes_enrichments -name "*fa" -type f | sort -h))
for genome in ${genomes[@]}; do

genome_base=${genome%.fa}
genome_base=${genome_base##*/}

# Clean up the name to remove any odd characters, to make more compatible downstream e.g., with MetAnnotate
genome_base=$(echo ${genome_base} | awk '{ gsub("[^A-Za-z0-9]", "_"); print }')

echo ${genome_base}
prokka --outdir ${work_dir}/prokka/${genome_base} --prefix ${genome_base} --locustag ${genome_base}_ --cpus 40 ${genome} 2&>/dev/null
# N.B., might have some problems with Archaea

done
```

## Part 2: ran MetAnnotate
Used Docker container version 0.9.2 relying on RefSeq from March 2017 (same as I've been using in the past)

### Started 181206 4:08AM JST. Finished in ~1-2 hours.
```
cd /Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/post-analysis/11_bin_env_abundance/vs2_dRep_and_enrichments
mkdir -p output/metannotate_contigs_vs1

enter-metannotate /Data/metannotate/data input/dereplicated_genomes/prokka input/hmm_files output/metannotate_contigs_vs1
ref_UID=$(stat -c "%u" /home/linuxbrew/output)
sudo chown -R linuxbrew /home/linuxbrew/output
echo 10 > MetAnnotate/concurrency.txt # set the number of threads
metannotate-wrapper-docker sequence orf_files hmm_files 2>&1 | tee output/metannotate_wrapper_docker.log
sudo chown -R $ref_UID /home/linuxbrew/output
# The 'chown' commands temporarily make the output folder belong to the linuxbrew user inside the Docker container so that the user can run the Docker commands. Files are given back to you at the end.
exit
```

# 7. Re-scanned raw reads with MetAnnotate to include the enrichment metagenomes and the new cyc2 HMMs.
First predicted ORFs for the metagenomic data using FGS++
Working in `/Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/post-analysis/11_bin_env_abundance/vs2_dRep_and_enrichments/input/metagenome_faa_predictions`

## Got hard links to existing predictions. Decided to just use R1 reads.
```
source_dir="/Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/post-analysis/03b_metannotate_raw_reads/input"
destination_dir="/Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/post-analysis/11_bin_env_abundance/vs2_dRep_and_enrichments/input/metagenome_faa_predictions"

files=($(find ${source_dir} -iname "*R1.frag.faa" | sort -h))

for file in ${files[@]}; do

ln ${file} ${destination_dir}/${file##*/}

done
```

## Installed FGS+ on mellea to closely mimic the install in neufeldserver (by running `conda list` there)
```
conda create -y -n FGS++ python=3.6.5
conda install -y -n FGS++ -c anaconda -c conda-forge meson=0.45.0 ninja=1.8.2 
conda install -y -n FGS++ libgcc-ng=7.2.0 libstdcxx-ng=7.2.0 certifi=2018.4.16 pip=10.0.1 setuptools=39.1.0 sqlite=3.23.1 wheel=0.31.1 openssl=1.0.2o
conda install -y -n FGS++ -c bioconda seqtk # this wasn't in the original environment but is needed

mkdir -p /home/jmtsuji/install
cd /home/jmtsuji/install
git clone https://github.com/unipept/FragGeneScanPlusPlus.git
# Commit 299cc18 - same as on neufeldserver

cd FragGeneScanPlusPlus

source activate FGS++

meson build
ninja -C build

# Then copy built binaries to conda env
cp build/FGS++ /home/jmtsuji/miniconda3/envs/FGS++/bin
cp -r train /home/jmtsuji/miniconda3/envs/FGS++/bin
# Then deleted install folder.
```

## Predicted genes of the R1 reads of the enrichments
Used FGS++ commit 299cc18 and seqtk version 1.3-r106
```
source activate FGS++

train_file_dir="/home/jmtsuji/miniconda3/envs/FGS++/bin/train"
threads=10
memory=10000

source_dir="/Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/post-analysis/11_bin_env_abundance/vs2_dRep_and_enrichments/input/metagenome_QC_reads/enrichments"
dest_dir="/Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/post-analysis/11_bin_env_abundance/vs2_dRep_and_enrichments/input/metagenome_faa_predictions"

input_files=(${source_dir}/L304_S_6D_QC_R1.fastq.gz ${source_dir}/L227_S_6D_QC_R1.fastq.gz)

for file in ${input_files[@]}; do
file_base=${file##*/}
file_base=${file_base%.fastq.gz}

echo ${file_base}

seqtk seq -A ${file} | FGS++ -s stdin -o ${dest_dir}/${file_base}.frag -w 0 -r ${train_file_dir} -t illumina_10 -p ${threads} -m ${memory}
done

conda list > ${dest_dir}/conda_env.log

```

## Ran MetAnnotate using Docker container version 0.9.2 relying on RefSeq from March 2017 (same as I've been using in the past)
Started 181206 at 5:17 JST
```
cd /Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/post-analysis/11_bin_env_abundance/vs2_dRep_and_enrichments
mkdir -p output/metannotate_raw_reads_vs1

enter-metannotate /Data/metannotate/data input/metagenome_faa_predictions input/hmm_files_small_set output/metannotate_raw_reads_vs1
ref_UID=$(stat -c "%u" /home/linuxbrew/output)
sudo chown -R linuxbrew /home/linuxbrew/output
echo 20 > MetAnnotate/concurrency.txt # set the number of threads
metannotate-wrapper-docker sequence orf_files hmm_files 2>&1 | tee output/metannotate_wrapper_docker.log
sudo chown -R $ref_UID /home/linuxbrew/output
# The 'chown' commands temporarily make the output folder belong to the linuxbrew user inside the Docker container so that the user can run the Docker commands. Files are given back to you at the end.
exit
```

# 8. Summarized output
Gathered output summary files together on my laptop at `~/Research_General/PhD/04b_Metagenome_resequencing_F2015/10_ATLAS_re_analysis/09_env_abundance_of_bins/vs2_dRep_scripted/output`

Had to combine a few files together because some genomes were added later. Was easy -- just bound the rows together (same columns).

Final summary files in `output/summary`

## Gathered assembly stats - percent assembled reads
# General command needed: `samtools view -c -F 4 ${bam_filename}`

```
bam_files=(/Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/coassembly_L227/CA-L227-2013/multi_mapping/L227-2013-6m.bam
/Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/coassembly_L227/CA-L227-2013/multi_mapping/L227-2013-8m.bam
/Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/coassembly_L227/CA-L227-2014/multi_mapping/L227-2014-6m.bam
/Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/coassembly_L227/CA-L227-2014/multi_mapping/L227-2014-8m.bam
/Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/coassembly_L442/CA-L442/multi_mapping/L442-2011-16-5m.bam
/Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/coassembly_L442/CA-L442/multi_mapping/L442-2014-15m.bam
/Analysis/jmtsuji/Documents/02_Active_workfolders/40_180721_GSB_genome_seq_analysis/L227_S_6D/sequence_alignment/L227_S_6D.bam
/Analysis/jmtsuji/Documents/02_Active_workfolders/40_180721_GSB_genome_seq_analysis/L304_S_6D/sequence_alignment/L304_S_6D.bam)

out_dir="/Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/post-analysis/11_bin_env_abundance/vs2_dRep_and_enrichments/output/percent_assembled_reads"
output_filepath="${out_dir}/assembled_read_stats.tsv"

# Initialize the output file
mkdir -p ${out_dir}
printf "metagenome\ttotal_reads\tmapped_assembled_reads\n" > ${output_filepath}

for bam_file in ${bam_files[@]}; do

filename_base=${bam_file%.bam}
filename_base=${filename_base##*/}

printf "${filename_base}\n"

printf "${filename_base}\t" >> ${output_filepath}
samtools view -c ${bam_file} >> ${output_filepath}
printf "\t" >> ${output_filepath}
samtools view -c -F 4 ${bam_file} >> ${output_filepath}
printf "\n" >> ${output_filepath}

done

```

# 9. Analyzed output
## Logic:
- Map to taxonomies (either via CheckM or via GTDB-tk). Add info to the table
- Further process stats, e.g., to average coverage normalized by total read count, or just mapped reads over total reads (or also normalize by genome length)
- Use MetAnnotate to find bins containing genes of interest. Filter to needed e-value cutoff. Those bins can then be showcased in their own panels.
- Make a table of the yes/no gene info for every bin. Map onto the main stats table.
- Plot as a bar plot:
- Overall: show bins > 1% of sample. Colour by similarity of taxonomic rank.
- dsrA: show ALL bins with dsrA
- cyc2: show ALL bins with cyc2


# LINK to file for mapping genomes between their GTDB and NCBI, etc. taxonomies: https://data.ace.uq.edu.au/public/gtdb/release86/bac_metadata_r86.tsv
# For release 86
# Downloaded to /Data/reference_databases/gtdbtk/release86 on mellea






## SIde note: aligned the cyc2 from other non-Chlorobi hits from metagenomes
### Manually grabbed from MetAnnotate contigs table
```
CA_L442_2__01967_hypothetical_protein_ma14_15_to_457
CA_L442_2__02700_hypothetical_protein_ma15_18_to_491
CA_L442_70__01747_hypothetical_protein_ma17_1_to_424
L227_2013_8m_001__03753_hypothetical_protein_ma23_25_to_478
L442_2011_16_5m_005__02494_hypothetical_protein_ma32_18_to_467
L227_2014_6m_013__00365_hypothetical_protein_ma24_13_to_456
L442_2011_16_5m_007__01101_hypothetical_protein_ma35_8_to_443
CA_L442_71__00809_hypothetical_protein_ma18_14_to_441
CA_L442_75__00831_hypothetical_protein_ma20_26_to_461
CA_L442_84__03042_hypothetical_protein_ma22_22_to_454
CA_L442_84__00357_hypothetical_protein_ma21_10_to_468

```

### Then parsed out the gene IDs and ran:
```
# Using
#clustalo 1.2.3
#seqtk 1.3-r106

ids=(CA_L442_2__01967
CA_L442_70__01747
L227_2013_8m_001__03753
L442_2011_16_5m_005__02494
L227_2014_6m_013__00365
L442_2011_16_5m_007__01101
CA_L442_71__00809
CA_L442_75__00831
CA_L442_84__03042)

#CA_L442_2__02700
#CA_L442_84__00357

work_dir="/Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/post-analysis/11_bin_env_abundance/vs2_dRep_and_enrichments/input/dereplicated_genomes/prokka"
id_file="${work_dir}/non_Chlorobia_cyc2.list"
out_file="${work_dir}/non_Chlorobia_cyc2_unaligned.faa"
aln_file_base="${work_dir}/non_Chlorobia_cyc2_aligned"

cd ${work_dir}
printf "" > ${out_file}

for id in ${ids[@]}; do
id_base=${id%__*}

echo ${id_base}
seqtk subseq ${id_base}/${id_base}.faa ${id_file} >> ${out_file}

done

clustalo -i ${out_file} --distmat-out=${aln_file_base}_distmat.txt --full --percent-id -o ${aln_file_base}.faa --threads 10 -v 2>&1 | tee ${aln_file_base}.log
```

