# QC metrics data Format

**Version: 1.0.0**

The qc metric data on the case and samples are recorded in a yaml outfile format with the following data structure and keys:


```
recipe: { #Hash of hashes
  [RECIPE_NAME]: {
    [METADATA_KEY] = string,
    version = string,
  },
},
sample: { #Hash of hashes
  [SAMPLE_ID]: {
    [INFILE_PREFIX]: {
      [PROGRAM_NAME]: {
        [HEADER_ID]: {
          [COMPARISON]: {
            [METRIC_KEY]: <data_type>,
          },
        },
        [METRIC_KEY]: <data_type>,
        version = string,
      },
    },
    [METRIC_KEY]: <data_type>,
  },
}
```
