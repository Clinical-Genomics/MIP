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
our $VERSION = 1.01;

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
        q{MIP::File::Path}     => [qw{ check_filesystem_objects_existance }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File::Path qw{ check_filesystem_objects_existance };

diag(   q{Test check_filesystem_objects_existance from Path.pm v}
      . $MIP::File::Path::VERSION
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

## Given a directory when it exists
my ($exist) = check_filesystem_objects_existance(
    {
        object_name    => $dir,
        object_type    => q{directory},
        parameter_name => q{dir},
    }
);

## Then return true
is( $exist, 1, q{Found directory} );

## Given a directory when it does NOT exist
($exist) = check_filesystem_objects_existance(
    {
        object_name    => q{does_not_exist},
        object_type    => q{directory},
        parameter_name => q{dir},
    }
);

## Then return false
is( $exist, 0, q{No directory} );

## Given a file when it exists
($exist) = check_filesystem_objects_existance(
    {
        object_name    => $file_name,
        object_type    => q{file},
        parameter_name => q{file_name},
    }
);

## Then return true
is( $exist, 1, q{Found file} );

## Given a file when it does not exist
($exist) = check_filesystem_objects_existance(
    {
        object_name    => q{does_not_exist},
        object_type    => q{file},
        parameter_name => q{file_name},
    }
);

## Then return false
is( $exist, 0, q{No file} );

done_testing();
