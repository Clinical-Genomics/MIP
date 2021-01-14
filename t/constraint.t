#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2018 };

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = ( q{MIP::Validate::Data} => [qw{ %constraint }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Validate::Data qw{ %constraint };

diag(   q{Test %constraint from Data.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a dir that does not exist
my $dir_does_not_exist = q{dir_does_not_exist};

## Then return false
is( $constraint{dir_exists}->($dir_does_not_exist),
    undef, q{Dir exists - directory does not exist} );

## Given a dir that exist

## Then return true
is( $constraint{dir_exists}->($Bin), 1, q{Dir exists - directory exists} );

## Given a file that does not exist
my $file_does_not_exist = q{file_does_not_exist};

## Then return false
is( $constraint{file_exists}->($file_does_not_exist), undef, q{File exists - file does not exist} );

## Given a file which exist
my $file_exist = catfile( $Bin, qw{ constraint.t } );

## Then return true
is( $constraint{file_exists}->($file_exist), 1, q{File exists - file exist} );

## Given a not executable file
my $file_is_not_executable = catfile( $Bin, qw{ constraint.t } );

## Then return false
is( $constraint{file_is_executable}->($file_is_not_executable),
    undef, q{File is executable - file is not executable} );

## Given an executable file
my $file_is_executable = catfile( $Bin, qw{ data modules miniconda envs mip_ci bin mip} );

## Then return true
is( $constraint{file_is_executable}->($file_is_executable),
    1, q{File is executable - file is executable} );

## Given a not gzipped file
my $filename = catfile(q{text.txt});

## Then return false
is( $constraint{is_gzipped}->($filename), 0, q{Is gzipped - plain file} );

## Given a gzipped file
my $gzipped_filename = catfile(q{test.gz});

## Then return true
is( $constraint{is_gzipped}->($gzipped_filename), 1, q{Is gzipped - gzipped file} );

## Given a dir that exist

## Then return false
is( $constraint{plain_file_exists}->($Bin), undef, q{Plain file exists - not a file} );

## Given a file that does not exist

## Then return false
is( $constraint{plain_file_exists}->($file_does_not_exist),
    undef, q{Plain file exists - file does not exist} );

## Given a file which exist

## Then return true
is( $constraint{plain_file_exists}->($file_exist), 1, q{Plain file exists - file exist} );

## Given a non-digit
my $non_digit = q{a_str};

## Then return false
is( $constraint{is_digit}->($non_digit), undef, q{Is digit - not a digit} );

## Given a digit

## Then return true
is( $constraint{is_digit}->(2), 1, q{Is digit - is a digit} );

done_testing();
