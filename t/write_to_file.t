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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_mip_hashes test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.01;

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
        q{MIP::Io::Write}      => [qw{ write_to_file }],
        q{MIP::Test::Fixtures} => [qw{ test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Io::Write qw{ write_to_file };

diag(   q{Test write_to_file from Write.pm v}
      . $MIP::Io::Write::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $test_dir = File::Temp->newdir();

## Given a hash when writing in YAML format
my %active_parameter = test_mip_hashes( { mip_hash_name => q{active_parameter}, } );
my $yaml_file_path   = catfile( $test_dir, q{write_yaml_to_file_test.yaml} );

write_to_file(
    {
        data_href => \%active_parameter,
        format    => q{yaml},
        path      => $yaml_file_path,
    }
);

## Then yaml file should exist
ok( -e $yaml_file_path, q{Created yaml file} );

## Given a hash when writing in TOML format
my $toml_file_path = catfile( $test_dir, q{write_toml_to_file_test.yaml} );

my %toml_hash = (
    hash  => q{toml},
    array => [qw{ toml toml }],
);

write_to_file(
    {
        data_href => \%toml_hash,
        format    => q{toml},
        path      => $toml_file_path,
    }
);

## Then yaml file should exist
ok( -e $toml_file_path, q{Created toml file} );

done_testing();
