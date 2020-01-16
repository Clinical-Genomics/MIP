#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir };
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
        q{MIP::Config}         => [qw{ set_default_config_dynamic_parameters }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Config qw{ set_default_config_dynamic_parameters };

diag(   q{Test set_default_config_dynamic_parameters from Config.pm v}
      . $MIP::Config::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given no supplied dynamic parameter
my @config_dynamic_parameters = qw{ analysis_constant_path };

my %active_parameter;

my %parameter = (
    analysis_constant_path => {
        default => q{analysis},
    }
);

set_default_config_dynamic_parameters(
    {
        parameter_href        => \%parameter,
        active_parameter_href => \%active_parameter,
        parameter_names_ref   => \@config_dynamic_parameters,
    }
);

## Then use default
is( $active_parameter{analysis_constant_path}, q{analysis}, q{Set default parameter} );

## Given a by user supplied dynamic parameter
$active_parameter{analysis_constant_path} = q{test_analysis};

set_default_config_dynamic_parameters(
    {
        parameter_href        => \%parameter,
        active_parameter_href => \%active_parameter,
        parameter_names_ref   => \@config_dynamic_parameters,
    }
);

is( $active_parameter{analysis_constant_path},
    q{test_analysis}, q{Set dynamic config parameter} );

done_testing();
