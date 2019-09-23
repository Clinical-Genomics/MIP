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
        q{MIP::Parse::Parameter} => [qw{ parse_dynamic_config_parameters }],
        q{MIP::Test::Fixtures}   => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Parse::Parameter qw{ parse_dynamic_config_parameters };

diag(   q{Test parse_dynamic_config_parameters from Parameter.pm v}
      . $MIP::Parse::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given

my %active_parameter = (
    cluster_constant_path  => catfile(qw{ root dir_1 dir_2 case_id! }),
    analysis_constant_path => q{analysis},
    case_id                => q{case_1},
    pedigree_file => catfile(qw{ cluster_constant_path! case_id!_pedigree.yaml }),
);

my @config_dynamic_parameters = qw{ cluster_constant_path analysis_constant_path };

my %parameter;

my $is_ok = parse_dynamic_config_parameters(
    {
        active_parameter_href         => \%active_parameter,
        config_dynamic_parameters_ref => \@config_dynamic_parameters,
        parameter_href                => \%parameter,
    }
);

## Then
ok( $is_ok, q{Parsed dynamic parameters and parmeters} );

done_testing();
