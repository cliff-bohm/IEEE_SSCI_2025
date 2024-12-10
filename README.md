This repository contains replication instructions for and additional results from "Towards Understanding Evolved Memory Retaining Circuitry in the Wireworld Cellular Automaton" presented at IEEE SSCI 2025

This work was run using the MABE (Modular Agent Based Evolver) framework. The two evolutionary envientments disscussed in the paper are included with the included MABE code (in code/Worlds).

# Additional Data
Additional data can be found in "ADDITIONAL_DATA". Image files as well as snapshot files that can be used to run analysis are avalible, organized by task (endlunds/nback) and replicate number (see "Replicating results" below for replication instructions).
Results are shown for all perfect performers generated during experiments.



# Replicating results

# MABE Installation Guide

## Step 1: Ensure Required Components (you do not need to install MABE yet)

Visit the [Installation Guide](https://github.com/Hintzelab/MABE/wiki/Installation-and-getting-started-with-MABE) to check required software:
- C++17 Compiler
- CMake >= 3.13.3
- Python >= 3.7

## Step 2: Build MABE

1. Download this repository.
2. Navigate to the top-level directory.
3. Run `sh tools/setup.cmd` to set up mbuild.
4. Run `./mbuild` (or `./mbuild.exe` for Windows).
5. The executable will be placed in the `work/` directory.

## Step 3: Run MABE

1. Navigate to the `work/` directory:
    ```sh
    cd work/
    ```
2. Run MABE using:
	./mabe (or ./mabe.exe) if on windows
	
# Running experiments

1. once you have a working MABE executable, cp the exactuable to either RUN_FILES/edlunds_5/. or RUN_FILES/nback_2/.
2. cd into the run directory you just placed MABE into
3. Execute MABE using:
    ```sh
    python ../tools/mq.py -l
    ```
	  mq.py is a wrapper tool that manages parameter settings.
      Parameters are set in the settings files (\*.cfg files) and mq_conditions.txt file (values in the mq_conditions file override).
      For further details on how `mq.py` works, refer to the [MQ documentation](https://github.com/Hintzelab/MABE/wiki/MQ).

      Information regaurding parameters can be found in `work/settings_world.cfg`
      To run different configurations changes can be made to `work/mq_conditions.txt`

4. when run, a directory will be created for the condition being run (for each task, we only used one condition).

	for edlunds maze this directory will be: C0__CBU_0__SX_10__SZ_10__NRV_0=
	for nback this directory will be:

	within this directory you will find a directory for each replicat run labeled with the random seed used for that replicat

The data files within each replicate data should be:

      max.csv: statistics for the individual in the population with the max score; recorded every 10 generations
      pop.csv: average statistics for all individual in the population; recorded every 10 generations
      snapshot_data_\*.csv: statistics for all organisms; recoded every 1000 generations
      snapshot_organisms_\*.csv: genomic data for all organisms; recoded every 1000 generations
	  
# processing information theory data
after completing a run


1. in order to generate visualizations, you need to switch over to the ANALYSIS directory
		for edlunds: cd ../../ANALYSIS/edlunds
		for nback:   cd ../../ANALYSIS/nback

 2. next you need to copy a snapshot_data_\*.csv and snapshot_organisms_\*.csv file to this directory
	these files should be renamed snapshot_data.csv and snapshot_organisms.csv

3. finally, run: sh run_commands.sh
      	this will first run mabe in visualize mode which will load the snapshot data and generate data related to the highest scoring agent in the population at the time of that snapshot. This make take a minute.
   	next the visualization code will be called.
   	you may need to install graphviz, details can be found here: https://graphviz.org/

4. cd to images/101 to see the final rendered images
5. you can change 101 to a different number in "run_commands.sh" in order to have final images sent to a different directory (the indention is to use the replicat number from the expiment).



## Useful Links

- [MABE GitHub Repository](https://github.com/Hintzelab/MABE/)
- [Installation and Getting Started Guide](https://github.com/Hintzelab/MABE/wiki/Installation-and-getting-started-with-MABE)
- [MQ Documentation](https://github.com/Hintzelab/MABE/wiki/MQ)
