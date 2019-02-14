# Sample information data format

**Version: 1.0.0**

The metadata on the case and samples are recorded in a yaml outfile format with the following data structure and keys:

```
analysis_date: string
case: string
[VCF_FILE_KEY]: {  #Hash of hashes
  clinical: {  
    path: string,
  },
  research: {  
    path: string,
  },
},
recipe: { #Hash of hashes
  [RECIPE_NAME]: {
    outdirectory: string,
    outfile: string,
    path: string,
    version: string,
    metafile_tag: {  
      directory: string,
      file: string,
      path: string,
      processed_by: string,
      version: string,
    },
  },
  [METAFILE_TAG]: {
    path: string,
  },
},
sample: { #Hash of hashes
  [SAMPLE_ID]: {
    recipe: {
      [RECIPE_NAME]: {
        outdirectory: string,
        outfile: string,
        path: string,
        version: string,
        metafile_tag: {
          directory: string,
          file: string,
          path: string,
          processed_by: string,
          version: string,
        },
      },
      [INFILE]: {
        [RECIPE_NAME]: {
          outdirectory: string,
          outfile: string,
          path: string,
          version: string,
          metafile_tag: {  
            directory: string,
            file: string,
            path: string,
            processed_by: string,
            version: string,
          },
        },
      },
      [METAFILE_TAG]: {  
        path: string,
      },
    },
  },
},
```
