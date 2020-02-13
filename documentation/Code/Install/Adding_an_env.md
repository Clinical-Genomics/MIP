# Adding a new environment

## Define the environment

Add the environment in the `templates/mip_install_<pipeline>_config_-<version>-.yaml` file.

1.  Add the environment name and define it

```Yaml
# [ENV_NAME] environment spec
e<env_name>:
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
      - <executable>
```

2.  Add the environment name to `definitions/install_<pipeline>_parameters.yaml`

```Yaml
e<env_name>:  
  associated_recipe:
    - mip
  data_type: HASH
  mandatory: no
  type: mip  
installations:
  associated_recipe:
    - mip
  data_type: ARRAY
  default:
    - e<env_name>
```

## Add the environment to the CLI

In `lib/MIP/Cli/Mip/Install/[PIPELINE]`:

1.  Add the new environment to the option `environment_name` under the `isa` key

2.  Add the new environment to the option `installations` under the `isa` key

## Add the environment to the travis tests

1.  Add the environment to the `load_env` option:

In `templates/mip_[PIPELINE]_config.yaml`:

```Yaml
load_env:
  mip_travis_<ENV_NAME>:
   <recipe | program>:
   installation: e<env_name>
   method: conda
```

2.  Add the executables of your environment to the test envs:

```Bash
$ mkdir -p t/data/modules/miniconda/envs/mip_travis_<env_name>
$ touch t/data/modules/miniconda/envs/mip_travis_<env_name>/bin/<program_executable>
$ chmod a+x t/data/modules/miniconda/envs/mip_travis_<env_name>/bin/<program_executable>
```
