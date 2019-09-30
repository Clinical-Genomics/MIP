Define parameters data format

**Version: 1.0.0**

The parameters that MIP support are recorded in a yaml file format. The order of parameters in this file does not matter for the execution of MIP.

**Rules**
- The definition file follows the yaml format.
- There are three types of parameters that can be defined: "recipe", "recipe_argument", "path".
- These keys are mandatory for each parameter: "type", "associated_recipe", "data_type".

The define parameters file has the following data structure and keys:
```
recipe: { # Recipe hash
  analysis_mode: string, value="sample | case" # Mode for running the analysis recipe
  associated_recipe: [
    - string, value="mip | <recipe>" # Evaluate parameter if associated recipe is active. Ties recipe argument to one or several recipes.
  ]
  data_type: string, value="SCALAR" # Parameter data type
  default: integer, value="0 (=off) | 1 (=on) | 2 (=simulate)"
  file_tag: string, value="nofile_tag | <file_tag>" # Tag to include in output filename or no tag use "nofile_tag"
  infile_suffix: string, # Suffix for infiles to analysis recipe
  outfile_suffix: string, # Set to enable file suffix for analysis recipe outfile
  program_executables: string, value="<executable_name>" # Set to check if can run when executing mip
  recipe_type: string, value="aligners | variant_callers | structural_variant_callers", # For collecting output from multiple analysis recipes
  type: string, value="recipe(=analysis_recipe_switch)"
}
recipe_argument: {
  associated_recipe: [
    - string, value="mip | <recipe>" # Evaluate parameter if associated recipe is active. Ties recipe argument to one or several recipes.
  ]
  data_type: string, value="SCALAR | ARRAY | HASH" # Parameter data type
  default: "SCALAR | ARRAY | HASH",
  element_separator: string, # Delimiter for input/output on cli
  is_reference: integer, value="1" # Defines this parameter as a reference to log
  mandatory: string, value="yes" # Do not supply to make optional
  type: string, value="mip(=global) | recipe_argument(=argument to a recipe)"
}
path_parameter: { # Path hash
  associated_recipe: [
    - string, value="mip | <recipe>" # Evaluate parameter if associated recipe is active. Ties recipe argument to one or several recipes.
  ]
  build_file: integer, value="0 | 1" # Build recipe switch, used to build reference meta data files
  data_type: string, value="SCALAR | ARRAY | HASH" # Parameter data type
  element_separator: string, # Delimiter for input/output on cli
  exists_check: integer, value="file | directory"  # Do not supply key to remove check
  is_reference: integer, value="1" # Defines this parameter as a reference to log
  mandatory: string, value="yes" # Do not supply to make optional
  reference: string, value="reference_dir", # Supply if parameter file object is expected in mip reference dir
  type: string, value="path(=file object)"
  update_path: string, value="absolute_path"
}
```
