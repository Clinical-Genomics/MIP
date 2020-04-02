#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Path qw{ rmtree };
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
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.04;

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
        q{MIP::Vcfanno}                    => [qw{ parse_toml_config_parameters }],
        q{MIP::Test::Fixtures}             => [qw{ test_log test_standard_cli }],
        q{MIP::Environment::Child_process} => [qw{ child_process }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Environment::Child_process qw{ child_process };
use MIP::Vcfanno qw{ parse_toml_config_parameters };

diag(   q{Test parse_toml_config_parameters from Vcfanno.pm v}
      . $MIP::Vcfanno::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Creates log object
my $log = test_log( {} );

## Replace file path depending on location - required for TRAVIS
my $test_reference_dir = catfile( $Bin, qw{ data references } );

### Prepare temporary file for testing
my $fqa_vcfanno_config =
  catfile( $test_reference_dir,
    qw{ grch37_frequency_vcfanno_filter_config_-v1.0-.toml  } );

# For the actual test
my $test_fqa_vcfanno_config = catfile( $test_reference_dir,
    qw{ grch37_frequency_vcfanno_filter_config_test_parse_toml_-v1.0-.toml  } );

my $file_path = catfile( $test_reference_dir, q{grch37_gnomad.genomes_-r2.0.1-.vcf.gz} );

## Replace line starting with "file=" with dynamic file path
my $parse_path =
    q?perl -nae 'chomp;if($_=~/file=/) {say STDOUT q{file="?
  . $file_path
  . q?"};} else {say STDOUT $_}' ?;

## Parse original file and create new config for test
my $command_string = join $SPACE,
  ( $parse_path, $fqa_vcfanno_config, q{>}, $test_fqa_vcfanno_config );

my %process_return = child_process(
    {
        commands_ref => [ $command_string, ],
        process_type => q{open3},
    }
);

## Given a toml config file
my %active_parameter = (
    frequency_filter   => 1,
    fqa_vcfanno_config => $test_fqa_vcfanno_config,
);

my $is_ok = parse_toml_config_parameters(
    {
        active_parameter_href => \%active_parameter,
    }
);

## Then return true
ok( $is_ok, q{Passed parsing for toml file} );

## Clean-up
rmtree($test_fqa_vcfanno_config);

done_testing();
