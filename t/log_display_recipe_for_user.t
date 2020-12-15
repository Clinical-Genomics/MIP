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
use MIP::Test::Fixtures qw{ test_log };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Log::MIP_log4perl} => [qw{ log_display_recipe_for_user }],
        q{MIP::Test::Fixtures}    => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Log::MIP_log4perl qw{ log_display_recipe_for_user };

diag(   q{Test log_display_recipe_for_user from MIP_log4perl.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given
my $indent_level = 2;
my $log          = test_log( {} );
my $recipe       = q{bwa_ar};

trap {
    log_display_recipe_for_user(
        {
            indent_level => $indent_level,
            log          => $log,
            recipe       => $recipe,
        }
    )
};

## Then
like( $trap->stderr, qr/INFO/xms, q{Throw info log message} );

done_testing();
