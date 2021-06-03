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
        q{MIP::List}           => [qw{ check_allowed_array_values }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::List qw{ check_allowed_array_values };

diag(   q{Test check_allowed_array_values from List.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given valid input
my @allowed_values = qw{ ok allowed permitted };
my @values         = qw{ ok permitted };

my $return = check_allowed_array_values(
    {
        allowed_values_ref => \@allowed_values,
        values_ref         => \@values,
    }
);

## Then return true
ok( $return, q{All elements are allowed} );

push @values, q{not_valid};

$return = check_allowed_array_values(
    {
        allowed_values_ref => \@allowed_values,
        values_ref         => \@values,
    }
);

## Then return false
is( $return, 0, q{Found not allowed element} );

done_testing();
