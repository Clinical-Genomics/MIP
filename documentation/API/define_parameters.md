Define parameters data format

**Version: 1.0.0**

The parameters that MIP support are recorded in a yaml file format. The order of parameters in this file does not matter for the execution of MIP.

**Rules**
- The definition file follows the yaml format.
- There are three types of parameters that can be defined: "pprogram", "program_argument", "program_path".
- These keys are mandatory for each parameter: "type", "associated_program", "data_type".

The define parameters file has the following data structure and keys:
```
pprogram: { # Program hash
  associated_program: [
    - string, value="mip | p<program>" # Evaluate parameter if associated program is active. Ties program argument to one or several programs.
  ]
  chain: string, values="MAIN | <chain> ", # Set the dependency tree chain for program parameter
  data_type: string, value="SCALAR" # Parameter data type
  default: integer, value="0 (=off) | 1 (=on) | 2 (=simulate)"
  file_tag: string, value="nofile_tag | <file_tag>" # Tag to include in output filename or no tag use "nofile_tag"
  infile_suffix: string, # Suffix for infiles to analysis recipe
  outdir_name: string, # Set outdirectory for program
  outfile_suffix: string, # Set to enable file suffix for analysis recipe outfile
  program_name_path: string, value="<binary_name>" # Set to check if can run when executing mip
  program_type: string, value="aligners | variant_callers | structural_variant_callers", # For collecting output from multiple analysis recipes
  type: string, value="program(=analysis_recipe_switch)"
}
program_argument: {
  associated_program: [
    - string, value="mip | p<program>" # Evaluate parameter if associated program is active. Ties program argument to one or several programs.
  ]
  data_type: string, value="SCALAR | ARRAY | HASH" # Parameter data type
  default: "SCALAR | ARRAY | HASH",
  element_separator: string, # Delimiter for input/output on cli
  mandatory: string, value="yes" # Do not supply to make optional
  type: string, value="mip(=global) | program_argument(=argument to a program)"
}
path_parameter: { # Path hash
  associated_program: [
    - string, value="mip | p<program>" # Evaluate parameter if associated program is activei. Ties program argument to one or several programs.
  ]
  build_file: integer, value="0 | 1" # Build recipe switch, used to build reference meta data files
  data_type: string, value="SCALAR | ARRAY | HASH" # Parameter data type
  element_separator: string, # Delimiter for input/output on cli
  exists_check: integer, value="file | directory"  # Do not supply key to remove check
  mandatory: string, value="yes" # Do not supply to make optional
  reference: string, value="reference_dir", # Supply if parameter file object is expected in mip reference dir
  type: string, value="path(=file object)"
  update_path: string, value="absolute_path"
}
```
