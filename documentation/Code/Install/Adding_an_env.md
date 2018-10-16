# Adding a new environment

## Define the environment
In the definitions folder add the environment to the install_[PIPELINE]_parameters.yaml.

1. Add the environment name and define it
"#" [ENV_NAME] environment spec
e[ENV_NAME]:
  conda:
    [PROGRAM_NAME]: [PROGRAM_VERSION]=[PROGRAM_VERSION_SUB_PATCH]
    pip:
    python: 3.6
  pip:
  shell:

2. Add the environment name to the `environment_name` hash
environment_name:
  e[ENV_NAME]

3. Add the environment name to the `installations` array
installations:
  - e[ENV_NAME]

## Add the environment to the CLI
In in `lib/MIP/Cli/Mip/Install/[PIPELINE]`

1. Add the new environment to the option `environment_name` under the `isa` key

2. Add the new environment to the option `installations` under the `isa` key
