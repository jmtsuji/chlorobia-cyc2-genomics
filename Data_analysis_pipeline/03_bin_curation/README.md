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
The output table here was cleaned up to produce `Table_S1.csv`, a copy of which is included in this folder.


## Manual cleanup of selected bins
All genome bins were imported into Anvi'o for subsequent manual cleanup of specific bins of interest. Bins of interest were seleted after the bin dereplication step based on taxonomy and CheckM stats (see manuscript).

### Lake metagenomes
Genome bins were imported into anvi'o 4 using [atlas-to-anvi.sh](https://github.com/jmtsuji/atlas-extensions), version ____.

Create the conda environment with all dependencies installed:
```
# TODO
```

Import the bins:
```
# User variables
input_dir="${github_repo_location}/Data_analysis_pipeline/03_bin_curation/01_dereplication"
output_dir="${github_repo_location}/Data_analysis_pipeline/03_bin_curation/02_anvio"
threads=12

# TODO
```

You can then refine bins of interest (e.g., `___`) by running the following code and then working with the anvi'o browser interface:
```
# TODO
```

Once finished, export the genome bin sequences using
```
# TODO
```

### Enrichment culture metagenomes
Because enrichment cultures were sequenced later, they were refined using updated software. Anvi'o 5 was used with [atlas-to-anvi.sh](https://github.com/jmtsuji/atlas-extensions), version ____.

Create the conda environment with all dependencies installed:
```
# TODO
```

Import the bins:
```
# User variables
input_dir="${github_repo_location}/Data_analysis_pipeline/03_bin_curation/01_dereplication"
output_dir="${github_repo_location}/Data_analysis_pipeline/03_bin_curation/02_anvio"
threads=12

# TODO
```

You can then refine bins of interest (e.g., `___`) by running the following code and then working with the anvi'o browser interface:
```
# TODO
```

Once finished, export the genome bin sequences using
```
# TODO
```


## Contig ordering via mauve
Contigs for the manually curated genomes were ordered based on the reference genome sequence of *Chlorobium luteolum*. The re-ordering does not actually affect the downstream analyses of this paper but is good practice in some cases (e.g., for comparing gene arrangement).

Installed the mauve contig mover (a bit difficult - not possible via conda):
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


