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
use Readonly;
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Io::Read qw{ read_from_file };
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
        q{MIP::Config}         => [qw{ check_cmd_config_vs_definition_file }],
        q{MIP::Io::Read}       => [qw{ read_from_file }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Config qw{ check_cmd_config_vs_definition_file };

diag(   q{Test check_cmd_config_vs_definition_file from Config.pm v}
      . $MIP::Config::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given no unique and hence illegal keys
my %parameter = read_from_file(
    {
        format => q{yaml},
        path   => catfile( $Bin, qw{ data test_data define_parameters.yaml } ),
    }
);

my %active_parameter = (
    bwa_mem                 => 1,
    vcfparser_outfile_count => 1,
    case_id                 => q{case_1},    #Add mandatory key default
    case_1                  => 1,
);

my $return = check_cmd_config_vs_definition_file(
    {
        active_parameter_href => \%active_parameter,
        parameter_href        => \%parameter,
    }
);

## Then return undef
is( $return, undef, q{No unique parameters} );

## Given illegal key
$active_parameter{illegal_key} = q{you shall not pass};

trap {
    check_cmd_config_vs_definition_file(
        {
            active_parameter_href => \%active_parameter,
            parameter_href        => \%parameter,
        }
    )
};

## Then fatal message should be thrown
like( $trap->stderr, qr/illegal\s+key/xms, q{Throw fatal message if illegal key} );

done_testing();
