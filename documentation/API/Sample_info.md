# Sample information data format

**Version: 1.0.0**

The metadata on the case and samples are recorded in a yaml outfile format with the following data structure and keys:

```
analysis_date: string
case: string
mip_version: string
[METAFILE_TAG]: {
  path: string,
},
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
},
sample: { #Hash of hashes
  [SAMPLE_ID]: {
    analysis_type: string
    [METAFILE_TAG]: {  
      path: string,
    },
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
        [INFILE]: {
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
    },
  },
},
```

## Methods
get_family_member_id:
Return hash with family member ids
```Perl
my %family_member_id = get_family_member_id(
    {
        sample_info_href => $sample_info_href,
    }
);
$family_member_id{children} = [<child1_id>, <child2_id>];
$family_member_id{father} = <father_id>;
$family_member_id{mother} = <mother_id>;
```

get_read_group:
Return hash with read group headers.
```Perl
my %read_group = get_read_group(
            {
                infile_prefix    => $infile_prefix,
                platform         => $active_parameter_href->{platform},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );
$rg{id} = <$infile_prefix>;
$rg{pu} = <flowcell>.<lane>.<sample_barcode>;
$rg{sm} = <$sample_id>;
$rg{pl} = <$platform>;
$rg{lb} = <$sample_id>; # Dummy value since the actual LB isn't available in MIP (yet)
```

get_sample_info_case_recipe_attributes:
Return case recipe attribute or attributes hash
```Perl
# Scalar
my $path = get_sample_info_case_recipe_attributes(
    {
        attribute        => q{path},
        recipe_name      => $recipe_name,
        sample_info_href => \%sample_info,
    }
);
$path = <string>;

# Hash
my %attribute = get_sample_info_case_recipe_attributes(
    {
        recipe_name      => $recipe_name,
        sample_info_href => \%sample_info,
    }
);
$attribute{path} = <string>;
```

get_sample_info_sample_recipe_attributes:
Return sample recipe attribute or attributes hash for infile key
```Perl
# Scalar
my $path = get_sample_info_sample_recipe_attributes(
    {
        attribute        => q{path},
        infile => $infile,
        recipe_name      => $recipe_name,
        sample_id => $sample_id,
        sample_info_href => \%sample_info,
    }
);
$path = <string>;

# Hash
my %attribute = get_sample_info_sample_recipe_attributes(
    {
        infile => $infile,
        recipe_name      => $recipe_name,
        sample_id => $sample_id,
        sample_info_href => \%sample_info,
    }
);
$attribute{path} = <string>;
```

get_sequence_run_type:
Return scalar sequence run type, (e.g. paired-end or single-end) or a hash of sequence run types per infile prefix
```Perl
# Scalar
my $sequence_run_type = get_sequence_run_type(
        {
            infile_lane_prefix => $infile_lane_prefix,
            sample_id          => $sample_id,
            sample_info_href   => $sample_info_href,
        }
    );
$sequence_run_type = <string>;

# Hash
my %sequence_run_type = get_sequence_run_type(
        {
            infile_lane_prefix_href => $infile_lane_prefix_href,
            sample_id               => $sample_id,
            sample_info_href        => $sample_info_href,
        }
    );
$sequence_run_type{$infile_lane_prefix} => <string>;
```

get_sequence_run_type_is_interleaved:
Return boolean depending on interleaved status of fastq file
```Perl
my $is_interleaved_fastq = get_sequence_run_type_is_interleaved(
        {
            infile_lane_prefix => $infile_prefix,
            sample_id               => $sample_id,
            sample_info_href        => $sample_info_href,
        }
    );
$is_interleaved_fastq = <boolean>;
```

set_file_path_to_store:
Set file path under store according to file type and file_tag
```Perl
set_file_path_to_store(
    {
        file_tag         => $file_tag,
        file_type        => $file_type,
        path             => $path,
        sample_info_href => \%sample_info,
    }
);
%sample_info = ( store => { $file_type => { $file_tag => $path, }, } );
```
