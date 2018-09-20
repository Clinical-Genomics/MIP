# Pedigree information data format

**Version: 1.0.0**

The metadata on the case is recorded in a yaml file format with the following data structure and keys:

```
family: string,
samples: [  # Array of hashes
  - { # Sample specific info hash
      analysis_type: string, value="wes | wgs | wts"
      capture_kit: string, 
      expected_coverage: Integer,
      father: string, value="<family_id> | 0"
      mother: string, value="<family_id> | 0"
      phenotype: string, value="affected | unaffected | unknown"
      sample_id: string,
      sex: string, value="male | female | other"
    }
  - {
      analysis_type: string, value="wes | wgs | wts"
    .,
    .,
    }
],
```
