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


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Script::Utils}  => [qw{ help }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Script::Utils qw{ help };

diag(   q{Test help from Utils.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given usage message
## Then print message to stdout and exit
trap {
    help(
        {
            exit_code => 1,
            usage     => q{testing},
        }
    )
};
is( $trap->stdout, qq{testing\n}, q{Print to STDOUT} );
ok( $trap->exit, q{Exit program} );

done_testing();
