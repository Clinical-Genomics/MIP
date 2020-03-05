# Tests
MIP uses the core perl module TEST::More, which is a framework for writing test.  All test scripts are located in the standard perl [t dir]. The entire test suite should be run prior to making a pull request. This is done by:
```Perl
prove t -r -j 4
```

## Test structure
When applicable the tests should follow the structure [Given-When-Then]. An example of test with such a structure can be found [here]

## Sub routines
Input parameters to each sub routine are validated via the core module validate params. Each sub routine should be accompanied by a test script checking the function of the sub routine.
**Template**: [test.t]

## Command modules
Each command module in MIP must be accompanied by a test for each sub routine that controls whether the command module yields the expected output.
**Template**: [command_module.t]

## Recipes
Are tested by test scripts and using sbatch traps upon execution. MIP thus relies on the processes coded within the recipe to supply a correct error code when throwing errors. This is catched through the traps and propagated to the STDERR file of the entire process or the sub process that threw the error. Since each command sub routine in MIP can write to the unix standard stream it is possible to add a stderr log file for each command process that is executed within the recipe.
**Template**: [recipe_analysis_case.t] [recipe_analysis_sample.t]

## Analysis
Are tested by running `mip_[PROCESS]_[PIPELINE].test` both locally and with continous integration using TRAVIS.
MIP will check that all files that are supposed to be produced by an analysis exist and have a file size larger than zero. Furthermore, several key qc parameters are evaluated and existence off certain vcf key value pairs are checked using Test::Harness and the mip_analysis.t script. All this is done in the analysis recipe analysisrunstatus.

## Templates
Basic templates for writing tests to new MIP modules is provided in [code dir] and have the standard `.t` file suffix.

## test.t
This is a basic template for writing tests for modules with sub routines that do not return a commands array and thus do not utilize the MIP/lib/MIP/Test/Commands.pm module for testing. The template includes a test to check the availability of the module that is to be tested.

## command_module.t
This template is intended for testing modules with sub routines that returns a `@commands` array, which consists of the base command followed by testable arguments. The arguments can be of three types either base, required or specific. The ´@command´ array is used by the test_function subroutine in the Commands module in order to control that the expected output is generated. Furthermore, the Commands module includes a test for writing to an already open filehandle.

Three types of input arguments can be tested using this template. Each type is hold in a separate hash array (`%base_argument`, `%required_argument` and `%specific_argument`). The required arguments are always supplied to the test function, either by themselves or in conjunction with a base or specific argument. The testable arguments can either be a Scalar ($scalar) or an Array (@array).

## recipe_analysis_case.t & recipe_analysis_sample.t
This is intended to test the code dependencies and compilation of the recipe and not the actual data output.

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
use MIP::Program::Conda qw{conda_install};
```  
#### Diagnose string
Change:
```Perl
diag(   q{Test SUB_ROUTINE from MODULE_NAME.pm v}
      . $MIP::PATH::TO::MODULE::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
. $EXECUTABLE_NAME );
```
To:
```Perl
diag(   q{Test conda_install from Conda.pm v}
      . $MIP::Program::Conda::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
. $EXECUTABLE_NAME );
```
#### Basic shell command
Change:
```Perl
my $function_base_command = q{BASE_COMMAND};
```
To:
```Perl
my $function_base_command = q{conda install};
```
#### Supplied arguments
##### Array input
Change:
```Perl
ARRAY => {
    inputs_ref      => [qw{ TEST_STRING_1 TEST_STRING_2 }],
    expected_output => q{PROGMRAM_OUTPUT},
},
```
To:
```Perl
packages_ref => {
    inputs_ref      => [qw{ test_package_1=1.2.3 test_package_2=1.2 }],
    expected_output => q{test_package_1=1.2.3 test_package_2=1.2},
},
```
##### Scalar input
Change:
```Perl
SCALAR => {
    input           => q{TEST_STRING},
    expected_output => q{PROGRAM_OUTPUT},
},
```
To:
```Perl
env_name => {
    input           => q{test_env},
    expected_output => q{--name test_env},
},
```
[t dir]: https://github.com/Clinical-Genomics/MIP/tree/develop/t
[Given-When-Then]: https://www.agilealliance.org/glossary/gwt/#q=~(filters~(postType~(~'page~'post~'aa_book~'aa_event_session~'aa_experience_report~'aa_glossary~'aa_research_paper~'aa_video)~tags~(~'given*20when*20then))~searchTerm~'~sort~false~sortDirection~'asc~page~1)
[here]: Tests/Exit_signals.md
[code dir]: https://github.com/Clinical-Genomics/MIP/tree/develop/templates/code/
[command_module.t]: https://github.com/Clinical-Genomics/MIP/tree/develop/templates/code/commad_module.t
[recipe_analysis_case.t]: https://github.com/Clinical-Genomics/MIP/tree/develop/templates/code/recipe_analysis_case.t
[recipe_analysis_sample.t]: https://github.com/Clinical-Genomics/MIP/tree/develop/templates/code/recipe_analysis_sample.t
[test.t]: https://github.com/Clinical-Genomics/MIP/tree/develop/templates/code/test.t
