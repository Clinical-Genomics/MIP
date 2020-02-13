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
use Path::Tiny qw{ path };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli };

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
        q{MIP::Test::Writefile} => [qw{ write_toml_config }],
        q{MIP::Test::Fixtures}  => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Test::Writefile qw{ write_toml_config };

diag(   q{Test write_toml_config from Writefile.pm v}
      . $MIP::Test::Writefile::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given test and toml paths
my $cluster_constant_path = catdir( dirname($Bin),          qw{ t data} );
my $test_reference_path   = catdir( $cluster_constant_path, q{references} );
my $toml_template_path    = catfile( $test_reference_path,
    q{grch37_frequency_vcfanno_filter_config_template-v1.0-.toml} );
my $toml_config_path =
  catfile( $test_reference_path, q{grch37_frequency_vcfanno_filter_config-v1.0-.toml} );

write_toml_config(
    {
        test_reference_path => $test_reference_path,
        toml_config_path    => $toml_config_path,
        toml_template_path  => $toml_template_path,
    }
);

## Then new config file should have been created
ok( -e $toml_config_path, q{Created test toml config path} );

## Create Path::Tiny object
$toml_config_path = path($toml_config_path);

## Read into memory
my $file_content = $toml_config_path->slurp_utf8;

my ($replaced_template_path) = $file_content =~ /$test_reference_path/ms;

## Then the created file should have replaced "TEST_REFERENCES!" in the config
ok( $replaced_template_path, q{Replaced template reference path} );

## Clean-up
unlink $toml_config_path;

done_testing();
