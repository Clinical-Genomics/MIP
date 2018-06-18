Initation file data format

**Version: 1.0.0**

The program that MIP supports per pipeline should be defined in the initiation file. It serves as a map to define the dependecies and order of execution of programs. This dependency tree is parsed to find the initaion point and dependencies of the `--start_with_program` option.

Programs (or analysis recipes switches) have dependencies for programs upstream within their chain. The `MAIN` chain is the general chain, which all other chains inherit from at different branching points. 

**Rules**:
- The initiation map follow the YAML-format.
- Keys in the definition map should have capital letters, except keys that are lowercase program names.
- Chains that inherit (or branch off) from the MAIN chain, is initiated with the key word `CHAIN_` followed by the ID of the chain, e.g. "CHAIN_PEDDY". These should be placed in the initation map at the point of separation according to the chains upstream dependency.
- Each program should be an element in an array.
- Each program name should be unique and only represented once in the initation map.

**Keys**
There are some keys that have special meaning in the initiation map:
```
ALL: Defines the complete dependency tree. Programs on this chain will have dependecies for all upstream programs independent of chain 
CHAIN_MAIN: Defines the main path through the analysis. This is the chain that all other chains should inherit from when branching the first time.
PARALLEL: Defines programs executed in parallel on seperate chains, but are then merged back into the origin chain after all parallel chains in the ´PARALLEL´ block has finished, e.g. in ensemble calling.
program_name: Within a `PARALLEL` block a key with lowercase matching a program name - defines a branching point run in parallel with more than one program within the new parallel chain that will be merged back to the origin chain once it has completed together with all other parallel chains. The first element within this parallel chain should match the lowercase program_name key. 
```

**Start_with_program**
- A program on the MAIN chain will have all downstream programs added to the list of programs to execute independent of chain.
- A program on a chain which is not the MAIN chain, will have all downstream programs within its own chain and the `ALL` chain added to the list of programs to execute.
- A program within a PARALLEL key will not add any other program from the PARALLEL list of program. Exceptions is when the program is within a lowercase program_name key. Then all downstream program within that list is added to the list of programs to execute. As usual downstream programs within the origin chain and `ALL` chain programs wil be added to the list to execute. 

**Reference**
rare_disease_initation.yaml(https://github.com/Clinical-Genomics/MIP/blob/develop/definitions/rare_disease_initiation.yaml)
rna_initiation.yaml(https://github.com/Clinical-Genomics/MIP/blob/develop/definitions/rna_initiation.yaml)


