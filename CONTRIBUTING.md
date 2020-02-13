## Versioning
MIP adheres to [semantic versioning](https://semver.org/).

## Branching model
MIP follows a gitflow [branching model](http://www.clinicalgenomics.se/development/dev/gitflow/).

## Code style
MIP uses Perl::Tidy, Perl::Critic and Yamllint in order to maintain readability and a consistent code style across the repo. These tools can be installed by running mip_install_perl.sh with the `-d` option when installing MIP's perl distribution and cpan modules.  
 
Prior to opening a pull request each commit should be tested with the supplied bash script `mip-check`. Running `bash mip-check` will test the code newly edited code using perltidy, perlcritic and yamllint.

 - Perltidy will automatically check and reformat your perl code. 
 - Perlcritic will lint the code and ensure that it follows perl best practice. 
 - Yamllint check any yaml files that have been altered. 

Try to keep the number of warnings to a minimum. This ensures that your new code adheres to the style that has been agreed upon for the MIP repo.
