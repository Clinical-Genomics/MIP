# Adding a PROGRAM

## Define the program and version
Add the program to a an env. In the definitions folder add the program to an environment in the install_[PIPELINE]_parameters.yaml.

- You need to add your program according to which method is used to install the program:

```
"#" [ENV_NAME] environment spec
e[ENV_NAME]:
  conda:
    [PROGRAM_NAME]: [PROGRAM_VERSION]=[PROGRAM_VERSION_SUB_PATCH]
  pip: [PROGRAM_NAME]=[PROGRAM_VERSION]
  shell: [PROGRAM_NAME]=[PROGRAM_VERSION]
```

## Add the environment to the CLI
In `lib/MIP/Cli/Mip/Install/[PIPELINE]`:

1. Add the new program to the option `select_programs` under the `isa` key

2. Add the new program to the option `skip_programs` under the `isa` key

3. Add the new program to the option `program_versions` under the `isa` key

4. Add the executable of your program to the test envs:

In `t/data/modules/miniconda/envs/`:
```
$ mkdir -p mip_travis_[ENV_NAME]
$ touch bin/[PROGRAM_EXECUTABLE]
$ chmod a+x bin/[PROGRAM_EXECUTABLE]
```
