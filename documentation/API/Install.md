# Install

**Version: 1.0**

## Singularity
```
singularity:
  [executable_name]:
    executable:
      [binary_1]: # If empty, expect binary to be in the PATH of the singularity image
      [binary_2]: /path/to/binary # Path to binary in image
      [placeholder]: "no_executable_in_image" # Use ENV variable in singularity image instead of explicit binary
    program_bind_paths: # Optional; Path to be bound to the container
      - /a/path
    uri: [singlularity_image_adress]
```
