# References

MIP can download many program prerequisites automatically via the mip download application ``mip download [PIPELINE]``. References available for download in MIP are located in: [definitions/define_download_references.yaml](https://github.com/Clinical-Genomics/MIP/blob/master/definitions/define_download_references.yaml) and the API is described [here](https://github.com/Clinical-Genomics/MIP/blob/develop/documentation/API/Download_references.md).

To add a new reference for download one needs to add a `reference_tag` and `version_tag` in the reference hash in the define_download_references.yaml file. Then you need to describe the meta data regarding your specific reference e.g the filename, url prefix, compression state etc, see the API for details.

To download you newly added reference using MIP, run:
```
$ mip download [PIPELINE] --reference_dir <reference_dir> --reference <reference_tag=version_tag>
```
This will create a batch file called `download_reference.sh`. To execute, run:
```Bash
$ bash download_reference.sh
```

To include the changes in the analysis - add your reference to your mip config.yaml file where appropriate.
