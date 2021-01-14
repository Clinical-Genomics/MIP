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
    my %perl_module = ( q{MIP::Active_parameter} => [qw{ update_dynamic_config_parameters }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Active_parameter qw{ update_dynamic_config_parameters };

diag(   q{Test update_dynamic_config_parameters from Active_parameter.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a dynamic case_id parameter i.e. "case_id!" in active_parameters
my %active_parameter = (
    analysis_constant_path => q{analysis},
    case_id                => q{case_1},
    cluster_constant_path  => catfile(qw{ root dir_1 dir_2 case_id! }),
);

## Given a cluster_constant_path when containing case_id!
my @dynamic_parameters = qw{ cluster_constant_path analysis_constant_path };

## When looping through all dynamic parameters and updating info
PARAMETER:
foreach my $parameter_name (@dynamic_parameters) {

    ## Updates the active parameters to particular user/cluster for dynamic config parameters following specifications. Leaves other entries untouched.
    update_dynamic_config_parameters(
        {
            active_parameter_ref   => \$active_parameter{$parameter_name},
            dynamic_parameter_href => { case_id => $active_parameter{case_id}, },
        }
    );
}

## Then cluster constant path should be updated with supplied case id
my $updated_cluster_constant_path = catdir(qw{ root dir_1 dir_2 case_1 });
my %expected_active_parameter     = (
    analysis_constant_path => q{analysis},
    case_id                => q{case_1},
    cluster_constant_path  => $updated_cluster_constant_path,
);
is_deeply( \%active_parameter, \%expected_active_parameter,
    q{Updated cluster constant path with case_id} );

done_testing();
