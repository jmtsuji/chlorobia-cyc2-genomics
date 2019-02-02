# ABOUT Figure 02 - Fe/S gene pathways among Chlorobia
Copyright Jackson M. Tsuji, Neufeld Research Group, 2019
Part of the larger *IISD-ELA Chlorobia cyc2 project*.

__TODO__: Gather the needed files into the riboprotein folder.

All code here is to be run in a Bash terminal of a unix system (e.g., Ubuntu 16 or 18) unless indicated otherwise.

# Choose your working directory -- you need to modify this variable based on where you download the Zenodo repository to:
For example:
```
WORK_DIR="/Analysis/jmtsuji/Chlorobia_cyc2_code/Figure_02_Fe_S_gene_pathways"
cd ${WORK_DIR}

# Also choose number of threads
THREADS=12
```

## 1. Made the phylogenetic tree
Conda environment for making the tree:
```
__TODO__ - check this works
conda create -y -n Chlorobia_cyc2_riboprotein_tree  -c conda-forge -c bioconda python=3.6 prokka=1.13.3 seqtk=1.3 clustalo=1.2.3 iqtree=1.6.7

# To run
source activate Chlorobia_cyc2_riboprotein_tree
```

### a. Gathered input genomes
#### Reference organisms
__TODO__

#### Genome bins from this study
__TODO__

### b. Ran a genome annotation on all genomes via prokka as a quick way to find the riboproteins
Made a table of the IDs of the riboproteins of interested, based on the work by Hug and colleagues (2016). See `riboprotein_tree/riboprotein_names.tsv`

```
mkdir -p ${WORK_DIR}/input_data/genomes/prokka
cd ${WORK_DIR}

genomes=($(find ${WORK_DIR}/input_data/genomes -maxdepth 1 -iname "*.fa" | sort -h))
for genome in ${genomes[@]}; do

genome_basename=${genome##*/}
genome_basename=${genome_basename%.*}
echo ${genome_basename}

prokka --outdir ${WORK_DIR}/input_data/genomes/prokka/${genome_basename} --prefix ${genome_basename} --cpus ${THREADS} ${genome}

done
```

### c. Did a quick automated scan over the genomes to see if any were lacking some riboprotein genes
```
output_dir=${WORK_DIR}/riboprotein_tree/riboproteins/search
mkdir -p ${output_dir}

curation=0
rp_names="${WORK_DIR}/riboprotein_tree/riboprotein_names.tsv"
rp_list="${WORK_DIR}/riboprotein_tree/riboprotein_names_prokka.list"
cut -d $'\t' -f 3 ${rp_names} | tail-n +2 > ${rp_list}

for genome in ${genomes[@]}; do

genome_basename=${genome##*/}
genome_basename=${genome_basename%.*}
echo ${genome_basename}

# Find all riboproteins in the prokka annotations file
cat ${rp_list} | xargs -I {} grep {}\$ ${WORK_DIR}/input_data/genomes/prokka/${genome_basename}/${genome_basename}.tsv | cut -d $'\t' -f 1,4,7 > ${output_dir}/${genome_basename}_riboproteins_IDs.tsv

# Get gene IDs for the riboproteins
cut -d $'\t' -f 1 ${output_dir}/${genome_basename}_riboproteins_IDs.tsv > ${output_dir}/${genome_basename}_riboproteins.list

# Iteratively grab the FAA entries for the genome, to ensure that they are kept in the intended order.
IDs=($(cat ${output_dir}/${genome_basename}_riboproteins.list))
printf "" > ${output_dir}/${genome_basename}_riboproteins.faa

for ID in ${IDs[@]}; do
echo ${ID} > ${output_dir}/${genome_basename}_tmp
seqtk subseq ${WORK_DIR}/input_data/genomes/prokka/${genome_basename}/${genome_basename}.faa ${output_dir}/${genome_basename}_tmp >> ${output_dir}/${genome_basename}_riboproteins.faa
done
${output_dir}/${genome_basename}_tmp

# Report a warning if not all were matched:
grep "^>" ${output_dir}/${genome_basename}_riboproteins.faa | cut -d ' ' -f 2- | sort -h > ${output_dir}/${genome_basename}_riboproteins_check.tmp

if [ $(head -n -1 ${rp_list} | sort -h | cmp - ${output_dir}/${genome_basename}_riboproteins_check.tmp > /dev/null; echo $?) = 1 ]; then
lines=$(cat ${output_dir}/${genome_basename}_riboproteins_IDs.tsv | wc -l)
echo "${name_base}: WARNING: Hit too few or too many riboproteins (${lines} total). Need to manually curate hits before proceeding."

# Move files to a curation_needed folder
mkdir -p ${output_dir}/curation_needed
mv ${output_dir}/${name_base}_riboproteins_IDs.tsv ${output_dir}/${name_base}_riboproteins.list ${output_dir}/${name_base}_riboproteins.faa ${output_dir}/curation_needed

# Set curation to 1 if at least 1 sample needs curation
curation=1

fi
rm ${output_dir}/${genome_basename}_riboproteins_check.tmp

done
```

For bins where curation is needed, there are two options:
- Duplicate genes were found -- here, one needs to be chosen.
- Missing genes -- here, either remove the genome from the list OR delete that gene from consideration in the gene tree

Looked through the entries manually and found that seven riboproteins were missing across all entries. Summarized those present across all genomes as `${WORK_DIR}/riboprotein_tree/riboproteins/shared_riboproteins.list`

__TODO__ did I also have to eliminate any duplicates??

Made list of riboproteins to ultimately use and saved as:
```
used_rp_file="${WORK_DIR}/riboprotein_tree/riboproteins/shared_riboproteins.list"
```

### d. Recovered shared riboproteins by repeating the above code but with ${used_rp_file} as the input list
```
output_dir=${WORK_DIR}/riboprotein_tree/riboproteins/shared_riboproteins
mkdir -p ${output_dir}

curation=0
rp_list="${WORK_DIR}/riboprotein_tree/riboproteins/shared_riboproteins.list"

for genome in ${genomes[@]}; do

genome_basename=${genome##*/}
genome_basename=${genome_basename%.*}
echo ${genome_basename}

# Find all riboproteins in the prokka annotations file
cat ${rp_list} | xargs -I {} grep {}\$ ${WORK_DIR}/input_data/genomes/prokka/${genome_basename}/${genome_basename}.tsv | cut -d $'\t' -f 1,4,7 > ${output_dir}/${genome_basename}_riboproteins_IDs.tsv

# Get gene IDs for the riboproteins
cut -d $'\t' -f 1 ${output_dir}/${genome_basename}_riboproteins_IDs.tsv > ${output_dir}/${genome_basename}_riboproteins.list

# Iteratively grab the FAA entries for the genome, to ensure that they are kept in the intended order.
IDs=($(cat ${output_dir}/${genome_basename}_riboproteins.list))
printf "" > ${output_dir}/${genome_basename}_riboproteins.faa

for ID in ${IDs[@]}; do
echo ${ID} > ${output_dir}/${genome_basename}_tmp
seqtk subseq ${WORK_DIR}/input_data/genomes/prokka/${genome_basename}/${genome_basename}.faa ${output_dir}/${genome_basename}_tmp >> ${output_dir}/${genome_basename}_riboproteins.faa
done
${output_dir}/${genome_basename}_tmp

# Report a warning if not all were matched:
grep "^>" ${output_dir}/${genome_basename}_riboproteins.faa | cut -d ' ' -f 2- | sort -h > ${output_dir}/${genome_basename}_riboproteins_check.tmp

if [ $(head -n -1 ${rp_list} | sort -h | cmp - ${output_dir}/${genome_basename}_riboproteins_check.tmp > /dev/null; echo $?) = 1 ]; then
lines=$(cat ${output_dir}/${genome_basename}_riboproteins_IDs.tsv | wc -l)
echo "${name_base}: WARNING: Hit too few or too many riboproteins (${lines} total). Need to manually curate hits before proceeding."

# Move files to a curation_needed folder
mkdir -p ${output_dir}/curation_needed
mv ${output_dir}/${name_base}_riboproteins_IDs.tsv ${output_dir}/${name_base}_riboproteins.list ${output_dir}/${name_base}_riboproteins.faa ${output_dir}/curation_needed

# Set curation to 1 if at least 1 sample needs curation
curation=1

fi
rm ${output_dir}/${genome_basename}_riboproteins_check.tmp

done
```

### e. Build the concatenated alignment
```
### Concatenate sequences (assume that FAA's are in order, as checked above)
source_dir="${WORK_DIR}/riboprotein_tree/riboproteins/shared_riboproteins"
output_dir="${WORK_DIR}/riboprotein_tree/riboproteins/alignment"
mkdir -p ${output_dir}

# Initialize output file
printf "" > ${output_dir}/Chlorobia_riboproteins_concatenated.faa

genome_basename=${genome##*/}
genome_basename=${genome_basename%.*}
echo ${genome_basename}

# Concatenate and rename. Also use 'sed' to get rid of any stop codons if present
(echo ">${genome_basename}" && grep -v "^>" ${source_dir_dir}/${genome_basename}_riboproteins.faa) | seqtk seq -A | sed -e 's/\*//g' >> ${output_dir}/Chlorobia_riboproteins_concatenated.faa

done

### Align
clustalo -i ${output_dir}/Chlorobia_riboproteins_concatenated.faa --seqtype=Protein --infmt=fa --outfmt=fa --threads=${THREADS} -v --force --wrap=100 -o ${output_dir}/Chlorobia_riboproteins_concatenated_aligned.faa > ${output_dir}/Chlorobia_riboproteins_concatenated_aligned_clustalo.log 2>&1
# Manually checked the alignment -- looks good.
```

### f. Build the tree
```
source_dir="${WORK_DIR}/riboprotein_tree/riboproteins/alignment"
output_dir="${WORK_DIR}/riboprotein_tree/riboproteins/phylogeny"
mkdir -p ${output_dir}

cd ${output_dir}
iqtree -s ${source_dir}/Chlorobia_riboproteins_concatenated_aligned.faa -o Ignavibacterium_album_JCM_16511 -pre Chlorobia_riboprotein_tree -nt AUTO -seed 63 -b 1000 -m MFP
# Takes several hours to run.
```
The tip names of the tree match the names of the genome files that were added.


## 2. Ran BackBLAST
### a. Collected predicted proteomes and query genes
For reference genomes, downloaded the predicted proteomes from NCBI:


For genome bins from this study, annotated using prokka previously (see ______)


Collected all sequences into `_____`. However, this folder is not saved in the repo to save space. You can repeat the code above yourself to gather the genomes if you'd like.
__TODO__

### b. Installed BackBLAST development commit cdaddc4
```
# First installed the dependencies via conda
conda create -y -n backblast_core -c bioconda python=2.7 biopython=1.72 blast=2.6.0

# Then got the script from the repo
wget https://github.com/LeeBergstrand/BackBLAST_Reciprocal_BLAST/archive/v1.0.tar.gz
tar -xvzf v1.0.tar.gz
rm v1.0.tar.gz
chmod 755 BackBLAST_Reciprocal_BLAST-1.0/BackBLAST.py

git clone https://github.com/LeeBergstrand/BackBLAST_Reciprocal_BLAST.git
cd BackBLAST_Reciprocal_BLAST
git checkout cdaddc4
cd ..

```

### c. Ran BackBLAST.py on all samples
```
source activate backblast_core

# Find relevant files
query_files=($(find ${WORK_DIR}/reference_genes -iname "*.faa" | sort -h))
ORF_files=($(find ${WORK_DIR}/input_ORFs -iname "*.faa" | sort -h))
e_value=1e-40
identity=20

mkdir -p backblast_results

# Run BackBLAST
echo "[ $(date -u) ]: Settings: e value cutoff = '${e_value}'; percent identity cutoff = ${identity}."
for query in ${query_files[@]}; do
query_basename=${query##*/}
query_basename=${query_basename%.faa}

# Get the name of the corresponding whole genome .faa file (HARD-CODED method based on file naming structure)
query_genome_basename=${query_basename%_gene_targets}
query_genome="${WORK_DIR}/input_ORFs/${query_genome_basename}.faa"

echo "[ $(date -u) ]: Searching for '${query_basename}' among the subjects. Using '${query_genome_basename}' as the reference predicted proteome."

for subject in ${ORF_files[@]}; do
subject_basename=${subject##*/}
subject_basename=${subject_basename%.faa}

output_filename="${WORK_DIR}/backblast_results/${query_genome_basename}__to__${subject_basename}.csv"
output_logname="${output_filename%.csv}.log"

echo "[ $(date -u) ]: 'BackBLAST.py -q ${query##*/} -r ${query_genome##*/} -s ${subject##*/} -e ${e_value} -i ${identity} -o ${output_filename##*/} 2>&1 | tee ${output_logname##*/}'"

${WORK_DIR}/BackBLAST_Reciprocal_BLAST/BackBLAST.py -q ${query} -r ${query_genome} -s ${subject} -e ${e_value} -i ${identity} -o ${output_filename} 2>&1 | tee ${output_logname}

done
done
echo "[ $(date -u) ]: Finished."
```
Copied screen output to a log at `backblast_loop.log`.

### d. Replaced C. ferro and C. clathratiforme with one-way BLAST due to issue with BackBLAST when BLAST'ing to self
When comparing a query to its own genome, BackBLAST currently seems to omit hits to the original gene -- maybe because of how the graph is constructed. This means that special measures are needed to compare a genome to itself. Should be pretty rudimentary; all hits should return 100% to self. No need to do reciprocal blast if the table looks normal.

#### Ran one-way BLAST on both queries to themselves
```
mkdir -p ${WORK_DIR}/backblast_results/one_way

query_file=${WORK_DIR}/reference_genes/Chl_ferrooxidans_KoFox_gene_targets.faa
ORF_file=${WORK_DIR}/input_ORFs/Chl_ferrooxidans_KoFox.faa
e_value=1e-40
identity=20 # Not actually needed here, but I kept this so I would remember to set the pident threshold later

blastp -query ${query_file} -subject ${ORF_file} -evalue ${e_value} -outfmt "10 qseqid sseqid pident evalue qcovhsp bitscore" > ${WORK_DIR}/backblast_results/one_way/Chl_ferrooxidans_KoFox__to__Chl_ferrooxidans_KoFox.csv

query_file=${WORK_DIR}/reference_genes/Chl_clathratiforme_BU_1_gene_targets.faa
ORF_file=${WORK_DIR}/input_ORFs/Chl_clathratiforme_BU_1.faa
e_value=1e-40
identity=20 # Not actually needed here, but I kept this so I would remember to set the pident threshold later

blastp -query ${query_file} -subject ${ORF_file} -evalue ${e_value} -outfmt "10 qseqid sseqid pident evalue qcovhsp bitscore" > ${WORK_DIR}/backblast_results/one_way/Chl_clathratiforme_BU_1__to__Chl_clathratiforme_BU_1.csv
```

#### Filtered using R to the best hits above the percent identity threshold (20)
See `backblast_results/one_way/filter_one_way_blast_hits.R` - run in interactive mode.

#### Replaced the original BackBLAST.py files with the new filtered one-way files
```
mkdir -p ${WORK_DIR}/backblast_results/old_reference
mv ${WORK_DIR}/backblast_results/Chl_ferrooxidans_KoFox__to__Chl_ferrooxidans_KoFox* ${WORK_DIR}/backblast_results/Chl_clathratiforme_BU_1__to__Chl_clathratiforme_BU_1* ${WORK_DIR}/backblast_results/old_reference
cp ${WORK_DIR}/backblast_results/one_way/Chl_ferrooxidans_KoFox__to__Chl_ferrooxidans_KoFox.csv.filtered ${WORK_DIR}/backblast_results/Chl_ferrooxidans_KoFox__to__Chl_ferrooxidans_KoFox.csv
cp ${WORK_DIR}/backblast_results/one_way/Chl_clathratiforme_BU_1__to__Chl_clathratiforme_BU_1.csv.filtered ${WORK_DIR}/backblast_results/Chl_clathratiforme_BU_1__to__Chl_clathratiforme_BU_1.csv
```
Now ready to proceed with plotting.

## 3. Plotted tree and heatmap using Figure_02_plotter.R
Required input files:
- `plot/gene_naming_info.tsv` - gives short abbreviations for all genes to be plotted and the plotting order. Also contains my personal notes.
- `plot/Chlorobia_naming_info.tsv` - links the tree tip names and backblast names to the final names to be plotted
- `riboprotein_tree/Chlorobia_riboprotein_tree.treefile` - the phylogenetic tree generated above
- `backblast_results/*.csv` -- the CSV files output from BackBLAST above.

Note: During iterative testing of the e-value cutoffs, I realized that a few genes had very poor hit profiles on reference organisms. They might be unreliable genes for reciprocal BLAST comparison. I commented them out of the final gene table so that they were excluded from the final plot:
- _qmoA_
- _cysA_
- _cysG_
- _dsrT_
- _soxJ_ and _soxK_ -- N.B., _soxJ_ is a cytochrome
- _dsrJ_
- Also the four proteins with unknown function but possible roles in H2 and sulfite reduction; these are not as relevant to the paper.

Ran `plot/Figure_02_plotter.R` to produce `plot/Figure_02_raw.pdf`. Note that you'll need to install all libraries loaded at the top of the script. After running, I then cleaned up the raw figure in Inkscape to make `plot/Figure_02_cleaned.pdf`, the final figure. Note that the script also outputs `plot/Figure_02_plotting_data.tsv` as a summary of the data to be plotted for the heatmap.

