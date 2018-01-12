# Modules
These comes in three flavors:
1. Perl core distribution modules that are installed when you install Perl
2. Cpanm perl modules usually installed using cpanm
3. MIPs internal perl modules that reside in within [lib]:

## Cpanm
All of the required cpanm perl modules in MIP are recorded in the standardised [cpanfile]:

This can be used to install or update cpanm modules using: 
```bash
cpanm --installdeps .
```
MIP will also read this file and check that all of these modules are installed prior to executing.

## MIPs internal modules
Standardised templates for writing MIPs internal modules e.g. command, generic module and recipe modules are found in the [code dir] with the standard perl `.pm` suffix.

## Using the templates
Change all occurrences of barewords that are written in capital letters to the appropriate variable. Notable exceptions are the words MIP, USAGE and VERSION.

### Examples
#### Paths to module
Change:
```Perl
use MIP::PATH::TO:MODULE qw{SUB_ROUTINE};
```
To:
```Perl
use MIP::Program::Qc::Fastqc;
```
#### Name of sub routine
Change:
```Perl
sub name_of_subroutine {
```
To:
```Perl
sub fastqc {
```
#### Basic shell command
Change:
```Perl
my $function_base_command = q{BASE_COMMAND};
```
To:
```Perl
my $function_base_command = q{fastqc};
```
[lib]: https://github.com/Clinical-Genomics/MIP/tree/master/lib/MIP/
[cpanfile]: https://github.com/Clinical-Genomics/MIP/tree/master/definitions/cpanfile
[code dir]: https://github.com/Clinical-Genomics/MIP/tree/master/templates/code/
