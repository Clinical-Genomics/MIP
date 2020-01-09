#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir };
use File::Temp qw{ tempdir tempfile };
use FindBin qw{ $Bin };
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
use MIP::Constants qw{ $COMMA $SPACE };
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
        q{MIP::Check::Path}    => [qw{ check_filesystem_objects_existance }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Check::Path qw{ check_filesystem_objects_existance };

diag(   q{Test check_filesystem_objects_existance from Path.pm v}
      . $MIP::Check::Path::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Create a temp dir
my $dir = tempdir( CLEANUP => 1 );

## Create a temp file using newly created temp dir
my ( $fh, $file_name ) = tempfile( DIR => $dir );

my %parameter = (
    dir       => $dir,
    file_name => { build_file => 1 },
);

### TEST
## Dirs
my ($exist) = check_filesystem_objects_existance(
    {
        parameter_name => q{dir},
        object_name    => $dir,
        object_type    => q{directory},
    }
);

is( $exist, 1, q{Found directory} );

($exist) = check_filesystem_objects_existance(
    {
        parameter_name => q{dir},
        object_name    => q{does_not_exist},
        object_type    => q{directory},
    }
);

is( $exist, 0, q{No directory} );

## Files
($exist) = check_filesystem_objects_existance(
    {
        parameter_name => q{file_name},
        object_name    => $file_name,
        object_type    => q{file},
    }
);

is( $exist, 1, q{Found file} );

($exist) = check_filesystem_objects_existance(
    {
        parameter_name => q{file_name},
        object_name    => q{does_not_exist},
        object_type    => q{file},
    }
);

is( $exist, 0, q{No file} );

done_testing();
