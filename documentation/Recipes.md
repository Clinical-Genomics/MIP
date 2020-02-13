# Recipes
Recipes are a list of ordered commands and subroutines that has a defined input and output. Many recipes can be shared across pipelines as they will inherit the input and set the output dynamically. For instance, pipeline `rd_dna` and `rd_dna_vcf_rerun` share many recipes. However, this is not always applicable depending on the pipeline run requirements (parallelization etc.).

## MIP currently has five recipe types:
1. Install
2. Download
3. Build
4. Analysis
5. Pipeline

Recipes templates for `install`, `download` and `analysis` can be found [here](https://github.com/Clinical-Genomics/MIP/blob/develop/templates/code):

The recipe parameter is also know as the analysis recipe switch. It turns on and off the recipe module.

To decide when and where the analysis recipe will be executed you add the name of the analysis recipe parameter to the initiation map file for the corresponding pipeline you are working on. The initiation map file is located in the `definition` folder in the MIP directory. The order of the initiation map and the chain you place the parameter on will decide when the recipe will execute and its upstream dependencies.

To add an analysis recipe or build recipe switch see [Parameters](https://github.com/Clinical-Genomics/MIP/blob/develop/documentation/Parameters.md):

Each recipes is attached to an environment. By default all recipes are reside in the MIP base environment, unless stated otherwise using the `load_env` parameter. This tells MIP to load a specific conda environment before executing any recipe specific shell instructions. MIP will check that all defined executables, (`program_executables` in define parameters), for that particular recipe and environment exists before launching submission to SLURM.
