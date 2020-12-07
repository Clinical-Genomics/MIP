## Versioning
MIP adheres to [semantic versioning].

## Branching model
MIP follows a gitflow [branching model].

## Dependencies
We use [Carton] to lock MIP's perl dependencies in the `cpanfile.snapshot`file.

### Initializing the environment
Carton will use the local directory to install modules into. You're recommended to exclude these directories from git by adding the `local/` to `.gitignore`. 
After installing Carton with `cpanm install Carton`, you can create/update the snapshot file by:

```
carton install
```

Once you've done installing all the dependencies, you can push your application directory to a remote machine (excluding local and .carton) and run the following command:

```
carton install --deployment
```

This will look at the cpanfile.snapshot and install the exact same versions of the dependencies into local, and now your application is ready to run.

The --deployment flag makes sure that carton will only install modules and versions available in your snapshot, and won't fallback to query for CPAN Meta DB for missing modules.

Another flavoe of this is [Carmel]. Unlike traditional CPAN module installer, Carmel keeps the build of your dependencies in a central repository, then select the library paths to include upon runtime.

## Code style
MIP uses Perl::Tidy, Perl::Critic and Yamllint in order to maintain readability and a consistent code style across the repo. These tools can be installed by running mip_install_perl.sh with the `-d` option when installing MIP's perl distribution and cpan modules.  
 
Prior to opening a pull request each commit should be tested with the supplied bash script `mip-check`. Running `bash mip-check` will test the code newly edited code using perltidy, perlcritic and yamllint.

 - Perltidy will automatically check and reformat your perl code. 
 - Perlcritic will lint the code and ensure that it follows perl best practice. 
 - Yamllint check any yaml files that have been altered. 

Try to keep the number of warnings to a minimum. This ensures that your new code adheres to the style that has been agreed upon for the MIP repo.

[branching model]: http://www.clinicalgenomics.se/development/dev/gitflow/
[Carmel]: https://metacpan.org/pod/Carmel
[Carton]: https://metacpan.org/pod/Carton
[semantic versioning]: https://semver.org/