#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Path qw{ remove_tree };
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
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Log::MIP_log4perl qw{ initiate_logger };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.03;

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
        q{MIP::Config}         => [qw{ write_mip_config }],
        q{MIP::Io::Read}       => [qw{ read_from_file }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Config qw{ write_mip_config };
use MIP::Io::Read qw{ read_from_file };

diag(   q{Test write_mip_config from Config.pm v}
      . $MIP::Config::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Creates log object
my $log = test_log( {} );

## Given a config file path
my %active_parameter = (
    config_file_analysis => catfile( $Bin, qw{ data test_data config_file } ),
    associated_program   => q{some_program},
    test_key             => q{This is a value},
);
my %sample_info;
trap {
    write_mip_config(
        {
            active_parameter_href => \%active_parameter,
            remove_keys_ref       => [qw{ associated_program }],
            sample_info_href      => \%sample_info,
        }
    )
};

## Then we should broadcast an INFO message
like( $trap->stderr, qr/INFO/xms, q{Send info log message} );

## Then some keys should be removed
is( $active_parameter{associated_program}, undef, q{Removed keys from active parameter} );

## Then the file should have been writen to disc and be reloaded
my %config_parameter = read_from_file(
    {
        format => q{yaml},
        path   => $active_parameter{config_file_analysis},
    }
);

is_deeply( \%config_parameter, \%active_parameter, q{Wrote config file} );

## Then the path of the config should be transferred to sample_info
is(
    $sample_info{config_file_analysis},
    $active_parameter{config_file_analysis},
    q{Transferred config file analysis path to sample_info}
);

## Clean-up
remove_tree( $active_parameter{config_file_analysis} );

done_testing();
