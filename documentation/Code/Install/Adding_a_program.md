# Adding a program

## Define the program and version
Add the program to an environment in the `templates/mip_install_<pipeline>_config_-<version>-.yaml` file.


- Add your program according to which method the installation process uses:

```Yaml
conda:
  <program_name>: <program_version>=<program_version_subpatch>
pip:
  <program_name>: <program_version>
shell:
  <program_name>:
    conda_dependency:
      <dependency>: <dependency_version>
    version: <program_version>
singularity:
  uri: <uri_to_program>
  executables:
    <executable>: <path to executable in container> | <blank if executable in container path>
```

## Add the program to the CLI
In `lib/MIP/Cli/Mip/Install/<pipeline>`:

1. Add the new program to the option `select_programs` under the `isa` key

2. Add the new program to the option `skip_programs` under the `isa` key

## Add program to integration test

  1. Add the executable of your program to the test envs:

```
$ touch t/data/modules/miniconda/envs/mip_ci/bin/<program_executable>
$ chmod a+x t/data/modules/miniconda/envs/mip_ci/bin/<program_executable>
```
2. Add the program to the template config

## Add installation tests
Add a path and/or a execution test to `templates/program_test_cmds.yaml`.
```Yaml
program_test_command:
  <program_name>:
    execution: '<program_command>'
    path: '<program_executable>'
```
The path test checks that the file is in your $PATH and is executable while the execution test executes a command and checks the exit status. Note that the command must return a zero exit status for the test to succeed.
