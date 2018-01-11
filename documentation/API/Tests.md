# Tests
Each command module in MIP should be accompanied by a test that controls whether the command module yields the expected output. MIP uses the cpan module TEST::More, which is a framework for writing test.
Basic templates for writing tests to new MIP modules is provided in [code dir] and have the standard `.t` file suffix:

## test.t
This is a basic template for writing tests for modules that do not return a commands array and thus do not utilize the MIP/lib/MIP/Test/Commands.pm module for testing. The template includes a test to check the availability of the module that is to be tested as well as the module required to display the help message.

## test_commands.t
This template is intended for testing modules that returns a `@commands` array, which consists of the base command followed by testable arguments. The arguments can be of three types either base, required or specific. The ´@command´ array is used by the test_function subroutine in the Commands module in order to control that the expected output is generated. Furthermore, the Commands module includes a test for writing to an already open filehandle.

Three types of input arguments can be tested using this template. Each type is hold in a separate hash array (`%base_argument`, `%required_argument` and `%specific_argument`). The required arguments are always supplied to the test function, either by themselves or in conjunction with a base or specific argument. The testable arguments can either be a Scalar ($scalar) or an Array (@array).

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
use MIP::Package_manager::Conda qw{conda_install};
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
      . $MIP::Package_manager::Conda::VERSION
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

[code dir]: https://github.com/Clinical-Genomics/MIP/tree/master/templates/code/
