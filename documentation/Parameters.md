# Parameters
All parameters in MIP are defined in the `definition` folder located in the MIP directory. The parameters are defined at the corresponding CLI command level YAML file:
 - "mip" i.e. "perl mip -> mip_parameters.yaml"
 - "[process]" e.g. "perl mip analyse -> analyse_parameters.yaml"
 - "[pipeline]" e.g. "perl mip analyse rare_disease -> rare_disease_parameters.yaml"

This makes MIP aware of the parameter and the keys and values used to define the parameter decides the features and methods that will be used on the parameter. You do not have to add the parameter in any specific order.

There are 3 types of parameters:
- "program"
- "program_argument"
- "path_parameter"

The program parameter is also know as the analysis recipe switch. It turns on and off the program module.

You should also have the parameter as an option on the command line interface. This is done by adding the parameter to the corresponding CLI perl module in `lib/MIP/Cli`.
 - "mip" i.e. "perl mip -> lib/MIP/Cli/Mip.pm"
 - "[process]" e.g. "perl mip analyse -> lib/MIP/Cli/Mip/Analyse.pm"
 - "[pipeline]" e.g. "perl mip analyse rare_disease -> lib/MIP/Cli/Mip/Analyse/Rare_disease.pm"

To decide when and where the analysis recipe will be executed you the name of the analysis recipe parameter to the initiation map file for the corresponding pipeline you are working on. The initiation map file is located in the `definition` folder in the MIP directory. The order of the initation map and the chain you place the parameter on will decide when the program will execute and its upstream dependencies.

So far we have defined the features of the analysis recipe, added it to the CLI, and decided where and when to execute the program module. We also have to point to the actual code of the recipe. This is done by adding your analysis recipe name and the code reference pointing to the analysis recipe perl module to the analysis_recipe and program_name hashes in the `lib/MIP/Recipe/Pipeline/[Pipeline]` perl module. 

To add a analysis recipe switch follow these steps:
 - Add the parameter at the corresponding CLI command level YAML file in the `definition` folder in the MIP directory
 - Add the parameter to the corresponing CLI perl module in `lib/MIP/Cli` in MIPs lib directory.
 - Add the analysis recipe switch to the initiation map file in the `definition` folder in the MIP directory.
 - Add the parameter name to the analysis_recipe and program_name hashes in the `lib/MIP/Recipe/Pipeline/[Pipeline]` perl module. 
