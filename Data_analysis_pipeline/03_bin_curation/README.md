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
- The Mauve Contig Mover must be installed manually (see instructions below)

## Bin dereplication
Unique genome bins were determined using dRep, version 2.0.5.  
dRep was run for the lake metagenomes. The enrichment culture genomes were determined later to be distinct from any of the lake metagenomes with an additional dRep run, but then the lake metagenome dRep results were used, with the enrichment culture genomes appended to the output table (because of the timing of when data became available for the manuscript).

### Create the conda environment with all needed dependencies:
```
conda create -y -n dRep -c bioconda -c r drep=2.0.5 r-dplyr
```
Use this environment via `conda activate dRep`

### Prepare CheckM and bin information for dRep
```
# User variables
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

### Dereplicate the bins
```
cd ${output_dir}/output
log_code=$(date '+%y%m%d_%H%M')
dRep dereplicate -p ${threads} -g ${output_dir}/input/bins/*.fasta --genomeInfo ${output_dir}/input/checkm_stats_reduced.csv ${output_dir}/output 2>&1 | tee dRep_${log_code}.log
```

Once dRep is done, make a summary file of dRep output and CheckM input using R  
First, navigate to the right folder and start R:
```
cd ${output_dir}
R
```

Then run the following R code:
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
The output table here was cleaned up to produce `Table_S1.csv`, a copy of which is included in the `01_dereplication` subfolder. Note that the enrichment bins had to be added manually.


## Manual cleanup of selected bins
All genome bins (pre-dereplication) were imported into Anvi'o for subsequent manual cleanup of specific bins of interest. Bins of interest were seleted after the bin dereplication step based on taxonomy and CheckM stats (see manuscript).

### Lake metagenomes
Genome bins were imported into anvi'o 4 using [atlas-to-anvi.sh](https://github.com/jmtsuji/atlas-extensions), version 1.0.22-coassembly-r4.

Create the conda environment with all dependencies installed:
```
conda create -y -n atlas_to_anvi_r4 -c bioconda -c conda-forge -c r anvio=4 diamond bwa bbmap gffutils r r-plyr r-dplyr r-getopt
git clone https://github.com/jmtsuji/atlas-extensions.git
cd atlas-extensions
git checkout 1.0.22-coassembly-r4
cd ..
cp atlas-extensions/parse*R atlas-extensions/atlas-to-anvi.sh ${CONDA_PREFIX}/bin
rm -rf atlas-extensions
```

The first time you use the environment, you'll need to download the COGs:
```
conda activate atlas_to_anvi_r4
anvi-setup-ncbi-cogs --just-do-it
```

Import the bins:
```
# User variables
atlas_dir="${github_repo_location}/Data_analysis_pipeline/02_assembly_and_binning/lake_metagenomes"
output_dir="${github_repo_location}/Data_analysis_pipeline/03_bin_curation/02_anvio"
sample_IDs=(CA-L227-2013 CA-L227-2014 CA-L442)
threads=12

echo "[ $(date -u) ]: importing ${#sample_IDs[@]} samples into anvi'o"

# Run samples in a loop:
for sample_ID in ${sample_IDs[@]}; do

    echo "[ $(date -u) ]: Importing sample '${sample}'"
    echo "[ $(date -u) ]: atlas-to-anvi.sh ${atlas_dir} ${sample_ID} ${output_dir}/${sample_ID} ${threads} 2>&1 | tee ${output_dir}/${sample_ID}.log"
    atlas-to-anvi.sh ${atlas_dir} ${sample_ID} ${output_dir}/${sample_ID} ${threads} 2>&1 | tee ${output_dir}/${sample_ID}.log

done

echo "[ $(date -u) ]: finished."
```

You can then refine bins of interest (e.g., `CA_L227_2013_55`) by running the following code and then working with the anvi'o browser interface:
```
bin_dir="${output_dir}/CA_L227_2013"
bin_name="CA_L227_2013_55"
taxonomic_level="t_genus"
cd ${bin_dir}
# Make a backup copy of the DBs to keep a record of pre-refinement bins
cp -r ${bin_dir##*/}_samples_merged ${bin_dir##*/}_samples_merged_BACKUP1_original
cp ${bin_dir##*/}_contigs.db ${bin_dir##*/}_contigs_BACKUP1_original.db
anvi-refine -p ${bin_dir##*/}_samples_merged/PROFILE.db -c ${bin_dir##*/}_contigs.db -C metabat2 -b ${bin_name} --taxonomic-level ${taxonomic_level} --title ${bin_name} -P 8080 --server-only

# This will stream the refinement interface to port 8080 on your server. You can then access the refinement interface on your PC's web browser by going to http://[server_name OR ip_address]:8080
# Alternatively, if running on your local machine, just remove '-P 8080 --server-only' and stream directly to your computer's web browser.
```

Once finished, export the genome bin sequences using
```
anvi-summarize -p ${bin_dir##*/}samples_merged/PROFILE.db -c ${bin_dir##*/}_contigs.db -C metabat2 -o ${bin_dir##*/}_summary_refined --taxonomic-level ${taxonomic_level} --init-gene-coverages 2>&1 | tee misc_logs/anvi-summarize-refined_lakes.log
```

I've included visual notes of how the bins were refined in `02_anvio/refinement_images`.
```

### Enrichment culture metagenomes
Because enrichment cultures were sequenced later, they were refined using updated software. Anvi'o 5 was used with [atlas-to-anvi.sh](https://github.com/jmtsuji/atlas-extensions), commit 99b85ac.

Create the conda environment with all dependencies installed:
```
conda create -y -n atlas_to_anvi_99b85ac -c bioconda -c conda-forge -c r anvio=5 diamond bwa bbmap gffutils r r-plyr r-dplyr r-getopt
git clone https://github.com/jmtsuji/atlas-extensions.git
cd atlas-extensions
git checkout 99b85ac
cd ..
cp atlas-extensions/parse*R atlas-extensions/atlas-to-anvi.sh ${CONDA_PREFIX}/bin
rm -rf atlas-extensions
```

The first time you use the environment, you'll need to download the COGs:
```
conda activate atlas_to_anvi_99b85ac
anvi-setup-ncbi-cogs --just-do-it
```

Import the bins:
```
# User variables
atlas_dir="${github_repo_location}/Data_analysis_pipeline/02_assembly_and_binning/enrichment_metagenomes"
output_dir="${github_repo_location}/Data_analysis_pipeline/03_bin_curation/02_anvio"
sample_IDs=(L227_S_6D L304_S_6D)
threads=12

echo "[ $(date -u) ]: importing ${#sample_IDs[@]} samples into anvi'o"

# Run samples in a loop:
for sample_ID in ${sample_IDs[@]}; do
    echo "[ $(date -u) ]: Importing sample '${sample}'"
    echo "[ $(date -u) ]: Creating mapping guide file"

    # Make the guide file
    mkdir -p "${output_dir}/mapping_guides"
    mapping_guide_filepath="${output_dir}/mapping_guides/${sample_ID}.tsv"
    raw_read_basepath="${atlas_dir}/sequence_quality_control/${sample_ID}_QC"
    printf "sample_name\traw_reads_filepaths\n" > ${mapping_guide_filepath}
    printf "${sample_ID}\t${raw_read_basepath}_R1.fastq.gz,${raw_read_basepath}_R2.fastq.gz,${raw_read_basepath}_se.fastq.gz\n" >> ${mapping_guide_filepath}

    echo "[ $(date -u) ]: Mapping guide file finished"

    echo "[ $(date -u) ]: Importing to anvio"
    echo "[ $(date -u) ]: atlas-to-anvi.sh normal ${atlas_dir} ${sample_ID} ${output_dir}/${sample_ID} ${threads} ${mapping_guide_filepath} 2>&1 | tee ${output_dir}/${sample_ID}.log"
    atlas-to-anvi.sh normal ${atlas_dir} ${sample_ID} ${output_dir}/${sample_ID} ${threads} ${mapping_guide_filepath} > ${output_dir}/${sample_ID}.log 2>&1
    echo "[ $(date -u) ]: Import finished"

done

echo "[ $(date -u) ]: Finished."
```

You can then refine bins of interest (e.g., `L227_S_6D_001`) by running the following code and then working with the anvi'o browser interface:
```
bin_dir="${output_dir}/L227_S_6D"
bin_name="L227_S_6D_001"
taxonomic_level="t_family"
cd ${bin_dir}
# Make a backup copy of the DBs to keep a record of pre-refinement bins
cp -r ${bin_dir##*/}_samples_merged ${bin_dir##*/}_samples_merged_BACKUP1_original
cp ${bin_dir##*/}_contigs.db ${bin_dir##*/}_contigs_BACKUP1_original.db
anvi-refine -p ${bin_dir##*/}_samples_merged/PROFILE.db -c ${bin_dir##*/}_contigs.db -C metabat2 -b ${bin_name} --taxonomic-level ${taxonomic_level} --title ${bin_name} -P 8080 --server-only

# This will stream the refinement interface to port 8080 on your server. You can then access the refinement interface on your PC's web browser by going to http://[server_name OR ip_address]:8080
# Alternatively, if running on your local machine, just remove '-P 8080 --server-only' and stream directly to your computer's web browser.
```

Once finished, export the genome bin sequences using
```
anvi-summarize -p ${bin_dir##*/}samples_merged/PROFILE.db -c ${bin_dir##*/}_contigs.db -C metabat2 -o ${bin_dir##*/}_summary_refined --taxonomic-level ${taxonomic_level} --init-gene-coverages 2>&1 | tee misc_logs/anvi-summarize-refined_enrichments.log
```

I've included visual notes of how the bins were refined in `02_anvio/refinement_images`.

## Contig ordering via mauve
Contigs for the manually curated genomes were ordered based on the reference genome sequence of *Chlorobium luteolum*. The re-ordering does not actually affect the downstream analyses of this paper but is good practice in some cases (e.g., for comparing gene arrangement).

Installed the mauve contig mover (a bit difficult - not possible directly via conda; **note that Java and seqtk are dependencies**):
```
work_dir="${github_repo_location}/Data_analysis_pipeline/03_bin_curation/03_contig_ordering"

mkdir -p "${work_dir}"
cd "${work_dir}"

## Downloaded the mauve java app
wget http://darlinglab.org/mauve/snapshots/2015/2015-02-13/linux-x64/mauve_linux_snapshot_2015-02-13.tar.gz
tar -xvzf mauve_linux_snapshot_2015-02-13.tar.gz
# This makes a folder called 'mauve_snapshot_2015-02-13'
rm mauve_linux_snapshot_2015-02-13.tar.gz
# Now Mauve.jar is at 'mauve_snapshot_2015-02-13/Mauve.jar'

## Downloaded progressiveMauve (also 2015-02-13 snapshot) and added to PATH
wget http://darlinglab.org/mauve/snapshots/2015/2015-02-13/linux-x64/progressiveMauve.bz2
bunzip2 progressiveMauve.bz2
chmod 755 progressiveMauve
sudo mv progressiveMauve /usr/local/bin # requires sudo!

# Prep for copying genomes
mkdir -p "${work_dir}/input_genomes"
```

Then, you must **manually** copy all refined genome bins into `${work_dir}/input_genomes`.

Download the *Chl. luteolum* reference genome
```
mkdir -p "${work_dir}/reference_genome"
cd "${work_dir}/reference_genome"

luteolum_fna_URL="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/012/485/GCF_000012485.1_ASM1248v1/GCF_000012485.1_ASM1248v1_genomic.fna.gz"
wget -O - ${luteolum_fna_URL} | gunzip > Chl_luteolum_DSM_273_genome.fa
```

Run the Mauve Contig Mover
```
# User variables
input_dir="${work_dir}/input_genomes"
output_dir="${work_dir}/ordered_genomes/original"
ref="${work_dir}/reference_genome/Chl_luteolum_DSM_273_genome.fa"
mauve_path="${work_dir}/mauve_snapshot_2015-02-13/Mauve.jar"

mkdir -p "${output_dir}" "${output_dir}../final"

draft_genomes=($(find ${input_dir} -maxdepth 1 -type f -iname "*.fa" | sort -h))

# Startup messages
echo "[ $(date -u) ]: Ordering contigs for ${#draft_genomes[@]} genomes"

for draft_genome in ${draft_genomes[@]}; do

# Set variables
base_name=${draft_genome%.*}
base_name=${base_name##*/}
output="${output_dir}/${base_name}"

# Run the contig ordering
mkdir -p ${output}
echo "[ $(date -u) ]: Running mauve contig mover on sample '${base_name}'"
java -Xmx1G -cp ${mauve_path} org.gel.mauve.contigs.ContigOrderer -output ${output} -ref ${ref} -draft ${draft_genome} > ${output}.log 2>&1

# Reformatted the output contigs (e.g., convert to uppercase)
ordered_filepath=$(find -mindepth 2 -maxdepth 2 "${output_dir}/${base_name}" -type f -name "${base_name}*.fas")
seqtk seq -A -U -l 100 "${ordered_filepath}" > "${output_dir}../final/base_name.fna"

done

echo "[ $(date -u) ]: Done."
```

## Final bin statistics
Calculated CheckM stats, length stats, and so on for the final set of genome bins

Installed dependencies:
```
conda create -n bin_stats -c bioconda bbmap checkm-genome prokka=1.13.3

# Must manually finish the checkm install by running once and specifying the database dir.
# Can use the 'checkm' folder in the ATLAS database folder from above
checkm
```
Activate this with `conda activate bin_stats` before running the code below.


```
out_dir="${github_repo_location}/Data_analysis_pipeline/03_bin_curation/04_bin_stats"
threads=12
final_bin_dir="${work_dir}/ordered_genomes/final"

mkdir -p ${out_dir}/rough ${out_dir}/prokka ${out_dir}/prodigal
cd ${out_dir}/rough

# Got length stats and such
statswrapper.sh in=${final_bin_dir}/*.fna > ../statswrapper.tsv

# Got checkM stats
checkm lineage_wf --file checkm/completeness.tsv --tab_table --quiet --extension fna --threads ${threads} ${final_bin_dir} checkm 2>&1 | tee checkm_lineage_wf.log
checkm tree_qa --tab_table --out_format 2 --file checkm/taxonomy.tsv checkm 2>&1 | tee checkm_tree_qa.log
cp checkm/completeness ../

# Annotated genomes and got tRNA counts
cd ${out_dir}/prokka
genome_bins=($(find ${final_bin_dir} -name "*.fna" | sort -h))

printf "Bin_ID\ttRNA_count\n" > ${out_dir}/tRNA_counts.tsv
for genome_bin in ${genome_bins[@]}; do

    bin_base=${genome_bin%.fna}
    bin_base=${bin_base##*/}

    prokka --outdir ${bin_base} --prefix ${bin_base} --cpus ${threads} ${genome_bin}

    num_tRNA=$(cut -d $'\t' -f 2 ${bin_base}/${bin_base}.tsv | grep "tRNA" | wc -l)
    printf "${bin_base}\t${num_tRNA}\n" >> ${out_dir}/tRNA_counts.tsv

done

# Predicted ORFs using prodigal at the same time, for reference (e.g., needed for Figure 2)
cd ${out_dir}/prodigal
for genome_bin in ${genome_bins[@]}; do

    bin_base=${genome_bin%.fna}
    bin_base=${bin_base##*/}

    prodigal -i ${genome_bin} -a ${bin_base}.faa -d ${bin_base}.ffn -f gff -o ${bin_base}.gff 2>&1 | tee ${bin_base}.log

done
```

The output tables were combined into the final Table 1 in the paper -- see Table 1 folder:
- statswrapper.tsv
- completeness.tsv
- tRNA_counts.tsv (note that the original counts were done manually, and the above coded method was made later for simplicity)

Now, the genome bins are done!


