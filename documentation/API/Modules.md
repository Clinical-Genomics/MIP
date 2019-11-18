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

### Command module
Each command module serves as an API for each specific command e.g. "cd", "git clone" or "bwa mem". The command module function is to standardise all interaction with the program in a versioned, testable way. Input parameters are validated through the validate params template in the each sub routine of the command module and output parameters are tested in the corresponding test script for each sub routine in the command module. The command module should be as simple and generic as possible i.e. contain no logic except to push or write commands depending on the input parameters. Modulation of the behavior of the sub routines is reserved for modulating the input parameters when calling the sub routine for instance from an analysis recipe.

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
use MIP::Program::Fastqc;
```
#### Name of sub routine
Change:
```Perl
our @EXPORT_OK = qw{ space separated subroutines };
.
.
sub name_of_subroutine {
```
To:
```Perl
our @EXPORT_OK = qw{ fastqc };
.
.
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
