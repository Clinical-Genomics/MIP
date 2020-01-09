#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use File::Temp;
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2018 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COLON $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Check::File}    => [qw{ check_mip_process_files }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Check::File qw{ check_mip_process_files };
use MIP::Unix::System qw{ system_cmd_call };

diag(   q{Test check_mip_process_files from File.pm v}
      . $MIP::Check::File::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $test_dir       = File::Temp->newdir();
my $test_file_path = catfile( $test_dir, q{test.sh} );

open my $filehandle, q{>}, $test_file_path
  or croak q{Cannot write to} . $SPACE . $test_file_path . $COLON . $SPACE . $OS_ERROR;

## Given path that exists and one that does not
my @paths = ( $test_file_path, catfile( cwd(), q{does_not_exist.test} ) );

my $is_ok = check_mip_process_files(
    {
        filehandle => $filehandle,
        paths_ref  => \@paths,
    }
);

close $filehandle;

## Then bash function should be written
ok( $is_ok, q{Wrote bash test} );

my %return = system_cmd_call( { command_string => q{bash } . $test_file_path } );

## Then "test.sh" file should be found
like( $return{output}[0], qr/Found\s+file/sxm, q{Found file test} );

## Then does_not_exist.test file should not be found
like( $return{error}[0], qr/Could\s+not\s+find/sxm, q{Could not find file test} );

done_testing();
