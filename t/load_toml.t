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
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Toml}           => [qw{ load_toml }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Toml qw{ load_toml };

diag(   q{Test load_toml from Toml.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a toml file to load
my $hash_ref = load_toml(
    {
        path => catfile(
            $Bin,
            qw{ data references grch37_frequency_vcfanno_filter_config_-v1.0-.toml }
        ),
    }
);

## Then hash should contain keys
ok( keys %{$hash_ref}, q{Loaded toml file} );

done_testing();
