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
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $DOLLAR_SIGN $FORWARD_SLASH $NEWLINE $SPACE };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

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
        q{MIP::File::Path}     => [qw{ check_allowed_temp_directory }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File::Path qw{ check_allowed_temp_directory };

diag(   q{Test check_allowed_temp_directory from Path.pm v}
      . $MIP::File::Path::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Creates log object
my $log = test_log( {} );

## Given not allowed temp dirs
my %active_parameter = ( outdata_dir => catfile(qw{a mip_analysis_outdata_dir}), );

my @is_not_allowed_temp_dirs = (
    $active_parameter{outdata_dir},
    $FORWARD_SLASH . q{scratch},
    $FORWARD_SLASH . q{scratch} . $FORWARD_SLASH,
);

TEMP_DIR:
foreach my $temp_dir (@is_not_allowed_temp_dirs) {

    trap {
        check_allowed_temp_directory(
            {
                not_allowed_paths_ref => \@is_not_allowed_temp_dirs,
                temp_directory        => $temp_dir,
            }
        )
    };

## Then exit and throw FATAL log message
    ok( $trap->exit, q{Exit if not allowed temp dir: } . $temp_dir );
    like( $trap->stderr, qr/FATAL/xms, q{Throw fatal log message: } . $temp_dir );

}

## Given allowed temp dirs
my $allowed_temp_dir =
  $FORWARD_SLASH . q{scratch} . $FORWARD_SLASH . $DOLLAR_SIGN . q{SLURM_JOB_ID};

my $is_allowed = check_allowed_temp_directory(
    {
        not_allowed_paths_ref => \@is_not_allowed_temp_dirs,
        temp_directory        => $allowed_temp_dir,
    }
);

## Then return "1"
is( $is_allowed, 1, q{Allowed dir: } . $allowed_temp_dir );

done_testing();
