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


## Constants
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

BEGIN {

    use Readonly;

## Constants
    Readonly my $COMMA   => q{,};
    Readonly my $NEWLINE => qq{\n};
    Readonly my $SPACE   => q{ };

    diag(   q{Test test_import from Fixtures.pm}
      . $COMMA
          . $SPACE . q{Perl}
          . $SPACE
          . $PERL_VERSION
          . $SPACE
          . $EXECUTABLE_NAME );

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Script::Utils}  => [qw{ help }],
        q{MIP::Test::Fixtures} => undef,
    );

    test_import( { perl_module_href => \%perl_module, } );
}

ok( 1, q{Passed imports} );

done_testing();
