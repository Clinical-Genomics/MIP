# MIP API

The most atomic unit in MIP is the `command modules` and internal sub routines. These modules are ordered in analysis recipes that in turn are assembled into pipeline recipes.

## Command Module
These modules serve as an API for each specific command e.g. "cd", "git clone" or "bwa mem". The command modules are independent, reusable, unit tested and ensure that all interaction with a specific command is recorded, versioned and tested in a single location in MIP. Every module can if supplied with an already opened filehandle write to a file or return an array of the commands depending on the input parameters. They can all write or append to the unix standard streams i.e. STDOUT, STDERR. They are always accompanied by a sub routine test script testing the expected input and output for each parameter to the command modules subroutines.

## Recipes
Recipes are a list of ordered commands and subroutines that has a defined input and output. They are always accompanied by a recipe test script testing the expected input and output and execution of sub routines within the recipe.

### MIP currently has five recipe types:
1. Install
2. Download
3. Build
4. Analysis
5. Pipeline

They are all located in `lib/MIP/Recipes/[PROCESS]`

### Install Recipes
Used in `lib/MIP/Main/Install.pm` to handle installation of programs and references.

### Download Recipes
Used in `lib/MIP/Main/Download.pm`to download and prepare references. 

### Build Recipes
Used in `lib/MIP/Main/Analyse.pm` to build meta datafiles for references.

### Analysis Recipes
Used in `lib/MIP/Main/Analyse.pm` to perform operation with a defined input and output for a process with specific upstream and/or downstream dependencies, e.g. qc of fastqc files or alignment and sorting of reads.

### Pipeline
Used in `lib/MIP/Recipes/[Pipeline].pm` to tie together several analysis recipes into entire pipelines.
