# Adding a new environment

## Define the environment
In the definitions folder add the environment to the install_[PIPELINE]_parameters.yaml.

1. Add the environment name and define it
```Yaml
# [ENV_NAME] environment spec
e[ENV_NAME]:
  conda:
    [PROGRAM_NAME]: [PROGRAM_VERSION]=[PROGRAM_VERSION_SUB_PATCH]
    pip:
    python: 3.6
  pip:
  shell:
```

2. Add the environment name to the `environment_name` hash
```Yaml
environment_name:
  e[ENV_NAME]
```

3. Add the environment name to the `installations` array
```Yaml
installations:
  - e[ENV_NAME]
```

## Add the environment to the CLI
In `lib/MIP/Cli/Mip/Install/[PIPELINE]`:

1. Add the new environment to the option `environment_name` under the `isa` key

2. Add the new environment to the option `installations` under the `isa` key

## Add the environment to the travis tests

1. Add the environment to the `load_env` option:

In `templates/mip_config.yaml`:
```Yaml
load_env:
  mip_travis_[ENV_NAME]:
   [RECIPE|PROGRAM]:
   installation: e[ENV_NAME]
   method: [ENV_METHOD] # e.g conda
```
2. Add the executables of your environment to the test envs:

In `t/data/modules/miniconda/envs/`:
```Bash
$ mkdir -p mip_travis_[ENV_NAME]
$ touch bin/[PROGRAM_EXECUTABLE]
$ chmod a+x bin/[PROGRAM_EXECUTABLE]
```
