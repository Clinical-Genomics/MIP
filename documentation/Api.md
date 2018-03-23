# MIP API

The most atomic unit in MIP is the command modules. These modules are ordered in analysis recipes that in turn are assembled into pipeline recipes.

## Command Module
These modules serve as an API for each specific command e.g. "cd", "git clone" or "bwa mem". The command modules are independent, reusable and ensure that all interaction with a specific command is recorded, versioned and tested in a single location in MIP. Every module can if supplied with an already opened filehandle write to a file or return an array of the commands depending on the input parameters. They can all write or append to the unix standard streams i.e. STDOUT, STDERR. They are always accompanied by a test script testing the expected input and output for each parameter to the command modules subroutines.

## Recipes
Recipes are a list of ordered commands and subroutines that has a defined input and output. MIP currently has 4 recipe types:
1. Install
2. Build
3. Analysis
4. Pipeline

### Install Recipes
Used in `Install.pm` to handle installation of programs and references.

### Build Recipes
Used in `Analyse.pm` to build metadatafiles for references.

### Analysis Recipes
Used in `Analyse.pm` to perform operation with a defined input and ouput for a process with specific upstream and/or downstream dependencies, e.g. qc of fastqc files or aligning and sorting reads.

### Pipeline
Used in `Analyse.pm` to string together several analysis recipes into entire pipelines.
