# Qc information data format

**Version: 1.0.0**

The QC data on the case and samples are recorded in a yaml outfile format with the following data structure and keys:

```
recipe: { #Hash of hashes
  [RECIPE_NAME]: {
    [KEY]: [VALUE],
    version: string,
  },
},
```

## Methods
set_qc_data_case_recipe_info:
Set case recipe arbitrary info in qc_data hash
```
set_qc_data_case_recipe_info(
    {
        key          => $key,
        qc_data_href => \%qc_data,
        recipe_name  => $recipe_name,
        value        => $value,
    }
);
my $value = $qc_data{recipe}{$recipe_name}{$key};
```
set_qc_data_case_recipe_version:
Set case recipe version in qc_data hash
```Perl
set_qc_data_case_recipe_version(
    {
        qc_data_href => \%qc_data,
        recipe_name  => $recipe_name,
        version      => q{a_version},
    }
);
```
