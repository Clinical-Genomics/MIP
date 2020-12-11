# Adding or updating  a program

1. Define the program and version in `templates/mip_install_config.yaml`

    ```Yaml
    container:
      <container>:
        bind_path:   #Optional
          <executable>: <path outside container>:<path inside container>
        executables:
          <executable>: <path to executable in container> | <blank if executable in container PATH> | "no_executable_in_image"
        uri: <uri_to_image>
    ```
    If you are only updating a program, skip to step 4

1. Add the program to the install CLI `lib/MIP/Cli/Mip/Install.pm`:

    1. Add the new program to the option `select_programs` under the `isa` key
    2. Add the new program to the option `skip_programs` under the `isa` key

1. Add the program to the correct pipeline (rd_dna/rd_rna) in`defintions/install_parameters`  

1. Install the program by running mip install --environment_name <mip_env> --select_program <your program>
