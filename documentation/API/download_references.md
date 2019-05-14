# References

**Version: 1.0.0**

The metadata on each reference is recorded in a yaml file format with the following data structure and keys:

```YAML
## Define reference tags
reference: # Hash of arrays
  reference_tag: [
    - version_tag,
  ]
reference_feature: { # Hash of hash
  reference_tag: { # Hash of hash
    genome_build_tag: { # Hash of hash
      version_tag: { # Hash
        file: <Name_of_file_@_origin>
        file_check: <Name_of_file_@_origin.check> # E.g. md5
        file_index: <Name_of_index_file_@_origin.index> # E.g. fai or csi
        file_index_check: <Name_of_index_file_@_origin> # E.g. md5
        outfile: <Name_of_outfile>
        outfile_check: <Name_of_outfile.check> # E.g. md5
        outfile_index: <Name_of_outfile.index> # E.g fai or csi
        outfile_index_check: <Name_of_outfile.check> # E.g. md5
        url_prefix: <url_prefix> # Up until filename
        outfile_decompress: decompress method # E.g. gzip
        outfile_index_decompress: decompress method # E.g. gzip
        outfile_check_method: check method # E.g md5sum
        outfile_index_check_method: check method # E.g md5sum
        outfile_reformat_command: reformat bash commands to apply to file post download
        outfile_bgzip_command: define bgzip command  # E.g to a downloaded decompressed file to index again
        outfile_tabix_command: define tabix command # E.g to index a downloaded decompressed, reformated and then compressed file again
      }
    }
  }
}
```
The complete file for the rd_dna pipeline can be found [here](https://github.com/Clinical-Genomics/MIP/blob/develop/templates/mip_download_rd_dna_config_-1.0-.yaml).
