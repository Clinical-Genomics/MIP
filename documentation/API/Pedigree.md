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
      father: string, value="<father_sample_id> | 0"
      subject_id: string, value="<subject_id>"
      mother: string, value="<mother_sample_id> | 0"
      phenotype: string, value="affected | unaffected | unknown"
      sample_id: string,
      sex: string, value="male | female | other"
      time_point: Integer,
    }
  - {
      analysis_type: string, value="wes | wgs | wts"
    .,
    .,
    }
],
```

#### Note
The key sample_id is unique for each molecular library. A subject can have many molecular libraries and the key subject_id in the pedigree file is used to link the sample_ids together. Some examples are shown below.
![Subject_id DNA-RNA][subject_id]  
<em>Here the subject_id connects samples across different analysis to one individual</em>  
![subject_id intra case][subject_id_2]  
<em>Subject_id can also connect samples from the same individual taken at different time points</em>

## Methods
get_pedigree_sample_id_attributes:
Return the value of for a supplied sample id with a given attribute (e.g. 'sex')
```Perl
my $sample_id_sex = get_pedigree_sample_id_attributes(
    {
        attribute => q{sex},
		sample_id => $sample_id,
		sample_info_href => $sample_info_href,
    }
);
```

[subject_id]: (Subject_id_description.png)
[subject_id_2]: (Subject_id_description_2.png)
