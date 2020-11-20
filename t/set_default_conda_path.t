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
use Test::Trap qw{ :stderr:output(systemsafe) };

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.02;

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
        q{MIP::Active_parameter} => [qw{ set_default_conda_path }],
        q{MIP::Test::Fixtures}   => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Active_parameter qw{ set_default_conda_path };

diag(   q{Test set_default_conda_path from Active_parameter.pm v}
      . $MIP::Active_parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

test_log( {} );

## Given a parameter named "conda_path"
my $parameter_name   = q{conda_path};
my %active_parameter = ( $parameter_name => undef, );

## When called with a faulty conda binary
trap {
    set_default_conda_path(
        {
            active_parameter_href => \%active_parameter,
            bin_file              => q{not_a_conda_path},
            conda_path            => $parameter_name,
        }
    );
};

## Then exit and throw FATAL log message
is( $trap->leaveby, q{exit}, q{Exit if the conda path cannot be found} );
like(
    $trap->stderr,
    qr/Failed \s+ to \s+ find \s+ path \s+ to /xms,
    q{Throw log message if no conda path}
);

## Given a proper call using "conda_path"

## When bin file is conda per default
set_default_conda_path(
    {
        active_parameter_href => \%active_parameter,
        conda_path            => $parameter_name,
    }
);

## Then a conda path should be set
ok( $active_parameter{$parameter_name}, q{Set conda path} );

done_testing();
