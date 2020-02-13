# Creating a container

This section describes how to create a singularity container using the template `templates/singularity/singularity_template.def`. Start by copying the template file, taking care to follow the naming convention.

```Bash
cp templates/singularity/singularity_template.def templates/singularity/Singularity.<program_name>-<program_version>
```

## Edit the file

The file has several headers which can be modified.

#### %help

Add a short help description of the container.

#### %labels

Add metadata such as maintainer and version of the definition file.

#### %environment

Here you define environment variables that should be available at runtime. Used similarly to how .bashrc is used.

#### %runscript

Add code that should execute when the container is invoked using either `singularity run` or `./<container_file>`.

#### %post

Code that is executed during the build process. Basic system programs can be installed by adding them to the list of programs to be installed via apt-get install. The template will download and set up conda so that you can install the programs needed via conda. You can also install via pip or build from source by adding to this section.
