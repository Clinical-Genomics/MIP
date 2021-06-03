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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Active_parameter} => [qw{ set_include_y }],

    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Active_parameter qw{ set_include_y };

diag(   q{Test set_include_y from Active_parameter.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a gender hash with males and others
my %active_parameter = (
    gender => {
        males   => [qw{ id_1 id_2 }],
        females => [qw{ id_3 }],
        others  => [qw{ id_4 }],
    },
);

## Then set include_y to 1
set_include_y(
    {
        active_parameter_href => \%active_parameter,
    }
);
is( $active_parameter{include_y}, 1, q{Set include_y to 1} );

## Given no males or other
delete $active_parameter{gender}{males};
delete $active_parameter{gender}{others};

## Then set include_y to 0
set_include_y(
    {
        active_parameter_href => \%active_parameter,
    }
);
is( $active_parameter{include_y}, 0, q{Set include_y to 0} );

done_testing();
