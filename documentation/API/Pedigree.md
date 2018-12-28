# Pedigree information data format

**Version: 1.0.0**

The metadata on the case is recorded in a yaml file format with the following data structure and keys:

```
case: string,
samples: [  # Array of hashes
  - { # Sample specific info hash
      analysis_type: string, value="wes | wgs | wts"
      capture_kit: string, 
      expected_coverage: Integer,
      father: string, value="<case_id> | 0"
      mother: string, value="<case_id> | 0"
      phenotype: string, value="affected | unaffected | unknown"
      sample_id: string,
      is_from_sample: string, value="<sample_id>""
      sex: string, value="male | female | other"
      time_point: Integer,
    }
  - {
      analysis_type: string, value="wes | wgs | wts"
    .,
    .,
    }
],
``

## Methods
get_pedigree_sample_id_attributes:
Return the value of for a supplied sample id with a given attribute (e.g. 'sex')
```Perl
my $sample_id_sex = get_pedigree_sample_id_attributes({ attribute => q{sex},
							sample_id => $sample_id,
							sample_info_href => $sample_info_href,
						      });
```
