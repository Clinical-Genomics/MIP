# Tests
Each command module in MIP should be accompanied by a test that controls whether the command module yields the expected output. MIP uses the cpan module TEST::More, which is a framework for writing test.
Basic templates for writing tests to new MIP modules is provided in:

    MIP/
    - templates/
        - code/
            - test.t
            - test_commands.t

## test.t
This is a basic template for writing tests for modules that do not return a commands array and thus do not utilize the MIP/lib/MIP/Test/Commands.pm module for testing. The template includes a test to check the availability of the module that is to be tested as well as the module required to display the help message.

## test_commands.t
This templates is for testing modules that returns a ´@commands´ array which, consists of the base command followed by testable arguments. The arguments can be of three types either base, required or specific. The ´@command´ array is used by the test_function subroutine in the Commands module in order to control that the expected output is generated. The required arguments are always supplied to the test function, either by themselves or in conjunction with a base or specific argument. The testable arguments can either be a Scalar or an Array.

## Using the templates
Change all occurrences of barewords that is written in capital letters to the appropriate variable. Notable exceptions are the words MIP, USAGE and VERSION.


### Examples
##### Paths to module
Change:
```Perl
use MIP::PATH::TO:MODULE qw{SUB_ROUTINE};
```
To:
```Perl
use MIP::PacketManager::Conda qw{conda_install};
```  
<br>
##### Basic shell command
Change:
```Perl
my $function_base_command = q{BASE_COMMAND};
```
To:
```Perl
my $function_base_command = q{conda install};
```
<br>
##### Supplied arguments
Change:
```Perl
my %required_argument = (
    ARRAY => {
        inputs_ref      => [qw{ TEST_STRING_1 TEST_STRING_2 }],
        expected_output => q{PROGMRAM_OUTPUT},
    },
);
```
To:
```Perl
my %required_argument = (
    packages_ref => {
        inputs_ref      => [qw{ test_package_1=1.2.3 test_package_2=1.2 }],
        expected_output => q{test_package_1=1.2.3 test_package_2=1.2},
    },
);
```
