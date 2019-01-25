# Initation map file data format

**Version: 1.0.0**

The recipe that MIP supports per pipeline should be defined in the initiation file. It serves as a map to define the dependencies, set chain IDs and order of execution of recipes. This dependency tree is parsed to find the initiation point and dependencies of the `--start_with_recipe` option.

Recipes (or analysis recipes switches) have dependencies for recipes upstream within their chain. The `MAIN` chain is the general chain, which all other chains inherit from at different branching points.

**Rules**:
- The initiation map follow the YAML-format.
- Keys in the definition map should have capital letters.
- Chains that inherit (or branch off) from the MAIN chain, is initiated with the key word `CHAIN_` followed by the ID of the chain, e.g. "CHAIN_PEDDY". These should be placed in the initiation map at the point of separation according to the chains upstream dependency. Exception to this rule is the dynamic generation of anonymous chain ID, which will have the recipe_name in capital letters and not include the key word `CHAIN_`.
- Each recipe should be an element in an array.
- Each recipe name should be unique and only represented once in the initiation map.

**Keys**
There are some keys that have special meaning in the initiation map:
```
CHAIN_ALL: Defines the complete dependency tree. Recipes on this chain will have dependencies for all upstream recipes independent of chain
CHAIN_MAIN: Defines the main path through the analysis. This is the chain that all other chains should inherit from when branching the first time.
PARALLEL: Defines recipes executed in parallel on separate chains, but are then merged back into the origin chain after all parallel chains in the ´PARALLEL´ block has finished, e.g. in ensemble calling.
recipe_name: Within a `PARALLEL` block a key with uppercase matching a recipe name - defines a branching point run in parallel with more than one recipe within the new parallel chain that will be merged back to the origin chain once it has completed together with all other parallel chains. The first element within this parallel chain should match the uppercase recipe_name key.
```

**Start_with_recipe**
- A recipe on the MAIN chain will have all downstream recipes added to the list of recipes to execute independent of chain.
- A recipe on a chain which is not the MAIN chain, will have all downstream recipes within its own chain and the `ALL` chain added to the list of recipes to execute.
- A recipe within a PARALLEL key will not add any other recipe from the PARALLEL list of recipe. Exceptions is when the recipe is within a uppercase recipe_name key. Then all downstream recipe within that list is added to the list of recipes to execute. As usual downstream recipes within the origin chain and `ALL` chain recipes will be added to the list to execute.

**Reference**
- [rare_disease_initation.yaml](https://github.com/Clinical-Genomics/MIP/blob/develop/definitions/rare_disease_initiation.yaml)
- [rna_initiation.yaml](https://github.com/Clinical-Genomics/MIP/blob/develop/definitions/rna_initiation.yaml)
