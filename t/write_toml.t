#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp;
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
use TOML::Tiny qw{ to_toml };

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COLON $COMMA $DOT $DOUBLE_QUOTE $NEWLINE $SINGLE_QUOTE $SPACE };
use MIP::Test::Fixtures qw{ test_mip_hashes test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

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
        q{MIP::Toml}           => [qw{ load_toml write_toml }],
        q{MIP::Test::Fixtures} => [qw{ test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Toml qw{ load_toml write_toml };

diag(   q{Test write_toml from Toml.pm v}
      . $MIP::Toml::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $test_dir       = File::Temp->newdir();
my $toml_file_path = catfile( $test_dir, q{ap_test.Toml} );

## Given a toml file to load
my $hash_ref = load_toml(
    {
        path => catfile(
            $Bin,
            qw{ data references grch37_frequency_vcfanno_filter_config_-v1.0-.toml }
        ),
    }
);

write_toml(
    {
        data_href => $hash_ref,
        path      => $toml_file_path,
    }
);

## Then Toml file should exist
ok( -e $toml_file_path, q{Created Toml file} );

## When loading toml file
my $written_toml_hash_ref = load_toml(
    {
        path => $toml_file_path,
    }
);

## Then hash from serialized toml should contain keys
ok( keys %{$written_toml_hash_ref}, q{Loaded serialized toml file} );

done_testing();
