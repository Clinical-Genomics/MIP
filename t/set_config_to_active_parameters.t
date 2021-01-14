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


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Config}         => [qw{ set_config_to_active_parameters }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Config qw{ set_config_to_active_parameters };

diag(   q{Test set_config_to_active_parameters from Config.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my %config = (
    hash   => { scalar => q{config_scalar_1}, },
    array  => [qw{ config_array_1 config_array_2 }],
    scalar => q{config_scalar_2},
);

my %active_parameter;

set_config_to_active_parameters(
    {
        config_parameter_href => \%config,
        active_parameter_href => \%active_parameter,
    }
);
## Tests

is( $active_parameter{hash}{scalar}, q{config_scalar_1}, q{Set hash from config} );

is( $active_parameter{array}[0], q{config_array_1}, q{Set array from config} );

is( $active_parameter{scalar}, q{config_scalar_2}, q{Set scalar from config} );

%active_parameter = (
    hash   => { scalar => q{active_scalar_1} },
    array  => [qw{ active_array_1 active_array_2 }],
    scalar => q{active_scalar_2},
);

set_config_to_active_parameters(
    {
        config_parameter_href => \%config,
        active_parameter_href => \%active_parameter,
    }
);

is( $active_parameter{hash}{scalar},
    q{active_scalar_1}, q{Did not set hash from config} );

is( $active_parameter{array}[0], q{active_array_1}, q{Did not set array from config} );

is( $active_parameter{scalar}, q{active_scalar_2}, q{Did not set scalar from config} );

done_testing();
