# References

MIP can download many program prerequisites automatically via the mip download application ``mip download [PIPELINE]``. Each pipeline has a yaml file where available references are specified ([rd_dna](https://github.com/Clinical-Genomics/MIP/blob/develop/templates/mip_download_rd_dna_config_-1.0-.yaml), [rd_rna]((https://github.com/Clinical-Genomics/MIP/blob/develop/templates/mip_download_rd_rna_config_-1.0-.yaml)) and the API is described [here](https://github.com/Clinical-Genomics/MIP/blob/develop/documentation/API/download_references.md).

To download the references using MIP, run:
```
$ mip download [PIPELINE] --reference_dir <reference_dir> --config <download_config>
```
This will submit SLURM jobs that automatically downloads all references specified in the config that is not already present in the reference directory.

For a full list of download options, run:
```bash
$ mip download [PIPELINE] --help
```
## Updating a reference
Open the download config corresponding to your pipeline and update or add a new `version_tag` in the reference hash. Update the same version tag together with any metadata changes in the reference_feature hash, see the API for details. Run the same mip command as above to download your updated reference.

To include the changes in the analysis - add your reference to your mip config.yaml file where appropriate.
