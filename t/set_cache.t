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
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli };

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
        q{MIP::Parameter}      => [qw{ set_cache }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Parameter qw{ set_cache };

diag(   q{Test set_cache from Set::Parameter.pm v}
      . $MIP::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a dynamic parameter
my %parameter = (
    test_dynamic => { update_path => q{absolute_path}, },
    test         => { not_dynamic => q{not_dynamic} },
);

set_cache(
    {
        parameter_href => \%parameter,
        aggregates_ref => [qw{ update_path:absolute_path secondary_key:does_not_exists }],
    }
);

my %expected_cache = ( absolute_path => [qw{test_dynamic}], );

## Then parameter should have been added to cash
is_deeply( \%{ $parameter{cache} }, \%expected_cache, q{Set cache parmeter to hash} );

## Given trailing garbage
trap {
    set_cache(
        {
            parameter_href => \%parameter,
            aggregates_ref => [q{update_path:absolute_path:trailing_garbage}],
        }
    )
};

## Then exit and throw FATAL log message
like(
    $trap->stderr,
    qr/Unexpected \s+ trailing \s+ garbage/xms,
    q{Throw fatal log message if trailing garbage is found}
);

done_testing();
