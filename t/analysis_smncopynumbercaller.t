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
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COLON $COMMA set_genome_build_constants $SPACE };
use MIP::Test::Fixtures qw{ test_log test_mip_hashes test_standard_cli };

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
        q{MIP::Recipes::Analysis::Smncopynumbercaller} =>
          [qw{ analysis_smncopynumbercaller }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Analysis::Smncopynumbercaller qw{ analysis_smncopynumbercaller };

diag(   q{Test analysis_smncopynumbercaller from Smncopynumbercaller.pm v}
      . $MIP::Recipes::Analysis::Smncopynumbercaller::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Constants
Readonly my $GENOME_VERSION_37 => 37;

my $log = test_log( { log_name => q{MIP}, no_screen => 1, } );

set_genome_build_constants(
    {
        genome_source  => q{grch},
        genome_version => $GENOME_VERSION_37,
    }
);

## Given analysis parameters
my $recipe_name    = q{smncopynumbercaller};
my $slurm_mock_cmd = catfile( $Bin, qw{ data modules slurm-mock.pl } );

my %active_parameter = test_mip_hashes(
    {
        mip_hash_name => q{active_parameter},
        recipe_name   => $recipe_name,
    }
);
$active_parameter{$recipe_name}                     = 1;
$active_parameter{recipe_core_number}{$recipe_name} = 1;
$active_parameter{recipe_time}{$recipe_name}        = 1;
my $sample_id = $active_parameter{sample_ids}[0];

my %file_info = test_mip_hashes(
    {
        mip_hash_name => q{file_info},
        recipe_name   => $recipe_name,
    }
);
%{ $file_info{io}{TEST}{$sample_id}{$recipe_name} } = test_mip_hashes(
    {
        mip_hash_name => q{io},
    }
);

my %infile_lane_prefix;
my %job_id;
my %parameter = test_mip_hashes(
    {
        mip_hash_name => q{recipe_parameter},
        recipe_name   => $recipe_name,
    }
);
@{ $parameter{cache}{order_recipes_ref} } = ($recipe_name);
my %sample_info;

my $is_ok = analysis_smncopynumbercaller(
    {
        active_parameter_href   => \%active_parameter,
        file_info_href          => \%file_info,
        infile_lane_prefix_href => \%infile_lane_prefix,
        job_id_href             => \%job_id,
        parameter_href          => \%parameter,
        profile_base_command    => $slurm_mock_cmd,
        recipe_name             => $recipe_name,
        sample_id               => $sample_id,
        sample_info_href        => \%sample_info,
    }
);

## Then return TRUE
ok( $is_ok, q{ Executed analysis recipe } . $recipe_name );

done_testing();
