# ABOUT assembly and binning
Copyright Jackson M. Tsuji, Neufeld Research Group, 2019
Part of the larger *IISD-ELA Chlorobia cyc2 project*.

All code here is to be run in a Bash terminal of a unix system (e.g., Ubuntu 16 or 18) unless indicated otherwise. Assumes you have already run code in `01_data_acquisition`. **WARNING: ONLY run this section of the code if you really want to perform several days' worth of sequencing assembly/annotation on your server!**

## Define where you downloaded the Zenodo repo:
```
zenodo_repo_location="/Analysis/jmtsuji/Chlorobia_cyc2_code"
```

## Depdendencies and requirements
- `docker` (e.g., `docker-ce`)
- An internet connection
- A server with at least ~300 GB storage and ~100 GB RAM!

## To set up the ATLAS coassembly docker container and associated wrappers (for ease of use) and databases
```
# Download the Docker container to use later
docker pull jmtsuji/atlas-extensions:1.0.22-coassembly-r3

# Download the Git repo and get the wrapper
cd ${zenodo_repo_location}/Data_analysis_pipeline/02_assembly_and_binning
wget https://github.com/jmtsuji/atlas-extensions/archive/1.0.22-coassembly-r3.tar.gz
tar -xzf 1.0.22-coassembly-r3.tar.gz
mv atlas-extensions-1.0.22-coassembly-r3/enter-atlas-coassembly .
rm -r atlas-extensions-1.0.22-coassembly-r3 1.0.22-coassembly-r3.tar.gz

# Download the databases
# **I recommend starting a screen/byobu/tmux session that you can pull out of after starting the download command -- will take several hours.
atlas_database_dir="${zenodo_repo_location}/Data_analysis_pipeline/02_assembly_and_binning/atlas_databases_1.0.22"
output_dir="${zenodo_repo_location}/Data_analysis_pipeline/02_assembly_and_binning/lake_metagenomes"
mkdir -p ${atlas_database_dir} ${output_dir}

cd ${zenodo_repo_location}/Data_analysis_pipeline/02_assembly_and_binning
./enter-atlas-coassembly ${atlas_database_dir} ${output_dir}
# You now enter the ATLAS container, like a virtual machine inside of your computer. You are root -- be careful!

atlas download --out-dir /home/atlas/databases
# This will take a while.

# When done, exit the docker container
exit

# The databases will now be saved at the ${atlas_database_dir} above.
```
Important note: sometime in 2017 or 2018, the conda install of prokka 1.12 stopped working. Prokka 1.13 works, however. If you run into a prokka error during ATLAS, you might need to modify the prokka version requirement number in `/opt/conda/envs/atlas_coassembly_env/lib/python3.6/site-packages/atlas/envs/required_packages.yaml` in the Docker container to `prokka=1.13`.

## Part A: Lake metagenomes
Get everything ready:
```
input_data_dir="${zenodo_repo_location}/Data_analysis_pipeline/01_data_acquisition/lake_metagenomes"
atlas_database_dir="${zenodo_repo_location}/Data_analysis_pipeline/02_assembly_and_binning/atlas_databases_1.0.22"
output_dir="${zenodo_repo_location}/Data_analysis_pipeline/02_assembly_and_binning/lake_metagenomes"

# Make the output directory for the assembly
mkdir -p ${output_dir}/tmp
cd ${output_dir}

# Note that a copy of the ATLAS config file used in the paper is present in the output directory: atlas_config_lake_metagenomes.yaml
# Please edit some of the settings if needed for your machine (e.g., threads, RAM). Note that the filepaths used for the databases and input files should already work if you've followed the code so far.
```

Perform the standard sequence assembly and binning workflow. **I recommend starting a screen/byobu/tmux session that you can pull out of after starting the download command -- could take several days.**
```
# Enter the Docker container. I have to ignore the typical wrapper command to include the input files:
docker run -v ${atlas_database_dir}:/home/atlas/databases \
        -v ${output_dir}:/home/atlas/output \
        -v ${input_data_dir}:/home/atlas/data \
        -it jmtsuji/atlas-extensions:1.0.22-coassembly-r3 /bin/bash
# You now enter the ATLAS container, like a virtual machine inside of your computer. You are root -- be careful!

# Start the run -- BUT please modify the number of jobss here to match the thread number you specified in the .yaml file, if you modified the .yaml file above!
config="/home/atlas/output/atlas_config_lake_metagenomes.yaml"

date_code=$(date +'%y%m%d_%H%M')
atlas assemble --dryrun --jobs 14 --out-dir /home/atlas/output ${config} 2>&1 | tee atlas_run_${date_code}_steps.log
atlas assemble --jobs 14 --out-dir /home/atlas/output ${config} 2>&1 | tee atlas_run_${date_code}.log
# Now wait a long time while the run finishes. Once done, exit the docker container:

exit
```

Now perform the three coassemblies, using standard assembly data as input. **I recommend starting a screen/byobu/tmux session that you can pull out of after starting the download command -- could take several days.** A wrapper around the ATLAS pipeline, [co-assembly.sh](https://github.com/jmtsuji/atlas-extensions/blob/master/co-assembly.sh) (version `1.0.22-coassembly-r3`), is used to perform co-assembly and is runnable via Docker container. The workflow is described in the Methods section of the paper. To start the Docker container:
```
atlas_database_dir="${zenodo_repo_location}/Data_analysis_pipeline/02_assembly_and_binning/atlas_databases_1.0.22"
output_dir="${zenodo_repo_location}/Data_analysis_pipeline/02_assembly_and_binning/lake_metagenomes"

cd ${output_dir}
./enter-atlas-coassembly ${atlas_database_dir} ${output_dir}
# You now enter the ATLAS container, like a virtual machine inside of your computer. You are root -- be careful!

# Configure the run. This generate a base config (.yaml) file for the run, but because the one from the paper will be used, you don't actually need this.
config=/home/atlas/output/atlas_config_lake_metagenomes_coassembly_TEMP.yaml
guide=/home/atlas/output/atlas_lake_metagenomes_coassembly_guide.tsv
co-assembly.sh /home/atlas/output /home/atlas/databases ${config} ${guide} 12

# Now perform the run 
config=/home/atlas/output/atlas_config_lake_metagenomes_coassembly.yaml
guide=/home/atlas/output/atlas_lake_metagenomes_coassembly_guide.tsv
log_code=$(date '+%y%m%d_%H%M')
co-assembly.sh /home/atlas/output /home/atlas/databases ${config} ${guide} 12 2>&1 | tee /home/atlas/output/co-assembly_${log_code}.log
# Now wait a long time while the run finishes. Output will be saved to /home/atlas/output/coassembly

# Once done, exit the docker container:
exit
```


## Part B: enrichment metagenomes
This feels a lot like working with the lake metagenomes above, except that only individual assemblies were used. Because he datasets are so much simpler, the assemblies should run much faster. Careful -- I worked with a different server this time, so the config file is set for 40 threads. Make sure your server has at least 40 threads, or else you'll need to modify the config.

Get everything ready:
```
input_data_dir="${zenodo_repo_location}/Data_analysis_pipeline/01_data_acquisition/enrichment_metagenomes"
atlas_database_dir="${zenodo_repo_location}/Data_analysis_pipeline/02_assembly_and_binning/atlas_databases_1.0.22"
output_dir="${zenodo_repo_location}/Data_analysis_pipeline/02_assembly_and_binning/enrichment_metagenomes"

# Make the output directory for the assembly
mkdir -p ${output_dir}/tmp
cd ${output_dir}

# Note that a copy of the ATLAS config file used in the paper is present in the output directory: atlas_config_enrichments.yaml
# Please edit some of the settings if needed for your machine (e.g., threads, RAM). Note that the filepaths used for the databases and input files should already work if you've followed the code so far.
```

Perform the standard sequence assembly and binning workflow. **I recommend starting a screen/byobu/tmux session that you can pull out of after starting the download command -- could take several hours.**
```
# Enter the Docker container. I have to ignore the typical wrapper command to include the input files:
docker run -v ${atlas_database_dir}:/home/atlas/databases \
        -v ${output_dir}:/home/atlas/output \
        -v ${input_data_dir}:/home/atlas/data \
        -it jmtsuji/atlas-extensions:1.0.22-coassembly-r3 /bin/bash
# You now enter the ATLAS container, like a virtual machine inside of your computer. You are root -- be careful!

# Start the run -- BUT please modify the number of jobss here to match the thread number you specified in the .yaml file, if you modified the .yaml file above!
config="/home/atlas/output/atlas_config_enrichments.yaml"

date_code=$(date +'%y%m%d_%H%M')
atlas assemble --dryrun --jobs 14 --out-dir /home/atlas/output ${config} 2>&1 | tee atlas_run_${date_code}_steps.log
atlas assemble --jobs 14 --out-dir /home/atlas/output ${config} 2>&1 | tee atlas_run_${date_code}.log
# Now wait a long time while the run finishes. Once done, exit the docker container:

exit
```

## Some final notes
- Your assemblies and bins will look different than in the paper because the algorithms are non-determininstic. Labels will also differ. Thus, the results are not directly cross-comparable.
- If you want the assemblies from the paper, they can be downloaded from NCBI: ___TODO___

Done! Next step is bin dereplication and cleanup.

