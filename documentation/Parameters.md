# Parameters
All parameters in MIP are defined in the `definition` folder located in the MIP directory. The parameters are defined at the corresponding CLI command level YAML file:
 - "mip" i.e. "perl mip -> mip_parameters.yaml"
 - "[process]" e.g. "perl mip analyse -> analyse_parameters.yaml"
 - "[pipeline]" e.g. "perl mip analyse rd_dna -> rd_dna_parameters.yaml"

This makes MIP aware of the parameter. The keys and values used to define the parameter decides the features and methods that will be used on the parameter. You do not have to add the parameter in any specific order.

There are 3 types of parameters:
- "recipe"
- "recipe_argument"
- "path"

The recipe parameter is also know as the analysis recipe switch. It turns on and off the recipe module. A feature of the analysis recipe switch is to decide if the recipe should analyse all samples or operate at the case level. The analysis mode is set as key "analysis_mode: [sample | case]" in the corresponding CLI command level YAML file in the definitions folder.

You should also have the parameter as an option on the command line interface. This is done by adding the parameter to the corresponding CLI perl module in `lib/MIP/Cli`.
 - "mip" i.e. "perl mip -> lib/MIP/Cli/Mip.pm"
 - "[process]" e.g. "perl mip analyse -> lib/MIP/Cli/Mip/Analyse.pm"
 - "[pipeline]" e.g. "perl mip analyse rd_dna -> lib/MIP/Cli/Mip/Analyse/Rd_dna.pm"

To decide when and where the analysis recipe will be executed you the name of the analysis recipe parameter to the initiation map file for the corresponding pipeline you are working on. The initiation map file is located in the `definition` folder in the MIP directory. The order of the initiation map and the chain you place the parameter on will decide when the recipe will execute and its upstream dependencies.

So far we have defined the features of the analysis recipe, added it to the CLI, and decided where and when to execute the recipe module. We also have to point to the actual code of the recipe. This is done by adding your analysis recipe name and the code reference pointing to the analysis recipe perl module to the analysis_recipe hash in the `lib/MIP/Recipe/Pipeline/[Pipeline]` perl module.

To add an analysis recipe switch follow these steps:
 - Add the parameter at the corresponding CLI command level YAML file in the `definition` folder in the MIP directory. Set the key `analysis_mode`to `sample` or `case`.
 - Add the parameter to the corresponding CLI perl module in `lib/MIP/Cli` in MIPs lib directory.
 - Add the analysis recipe switch to the initiation map file in the `definition` folder in the MIP directory.
 - Add the parameter name to the `analysis_recipe` hash in the `lib/MIP/Recipes/Pipeline/[Pipeline]` perl module.
 - Fill your analysis_recipe with content.

To add a build recipe switch follow these steps:
 - Add the parameter at the corresponding CLI command level YAML file in the `definition` folder in the MIP directory. Set the keys `build_file` to `1`, `mandatory` to `no` and connect it to one or more recipes using the key `associated_recipe`.
 - Add the parameter to the file_info hash in the corresponding CLI perl module in `lib/MIP/Cli` in MIPs lib directory (add input parameters used in the build recipe to the CLI if required).
 - If custom defaults are required like setting the human genome reference per default depending on the `human_genome_reference` option. Add the build parameter in `lib/MIP/MAIN/Analyse.pm` for the call to sub `set_custom_default_to_active_parameter`. Add the build parameter to the `set_to_active_parameter` hash in `MIP/Set/Parameter` for sub `set_custom_default_to_active_parameter` and add the `_set_human_genome` code reference.
 - Add the build_recipe coderef to the build_recipe hash in the corresponding pipeline in `lib/MIP/Recipes/Build/[PIPELINE]` perl module.
 - Fill your build_recipe with content.

To add a download recipe switch follow these steps:
 - Add the parameter at the corresponding CLI command level YAML file in the `definition` folder in the MIP directory. Set the key `analysis_mode`to `case`.
 - Add the reference and version tags to the reference hash in `template/mip_download_[Pipeline]_-[Version]-.yaml`. In the same file define what to download and how by adding the reference, genome_build, reference tags and then keys according to the [download API](https://github.com/Clinical-Genomics/MIP/blob/develop/documentation/API/download_references.md).
 - Add the same reference info to the test fixture `t/data/test_data/download_active_parameters.yaml` for your download_[recipe_tag].t script.
 - Add the download_recipe coderef to the `download_recipe` hash in the corresponding pipeline `lib/MIP/Recipes/Pipeline/[Pipeline]` perl module.
 - Fill your download_recipe with content.
