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
use MIP::Test::Fixtures qw{ test_log test_mip_hashes test_standard_cli };

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
        q{MIP::Analysis}       => [qw{ update_recipe_mode_for_pedigree }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Analysis qw{ update_recipe_mode_for_pedigree };

diag(   q{Test update_recipe_mode_for_pedigree from Analysis.pm v}
      . $MIP::Analysis::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { no_screen => 1, } );

## Given pedigree and recipe
my %sample_info = test_mip_hashes( { mip_hash_name => q{qc_sample_info}, } );
my %active_parameter = ( blobfish => 1 );

## When pedigree phenotypea are compatible
update_recipe_mode_for_pedigree(
    {
        active_parameter_href => \%active_parameter,
        recipes_ref           => [qw{ blobfish }],
        sample_info_href      => \%sample_info,
    }
);

## Then keep blobfish recipe active
is( $active_parameter{blobfish}, 1, q{Recipe is compatible with pedigree} );

## When non compatible pedigree phenotypes
$sample_info{sample}{ADM1059A1}{phenotype} = q{unaffected};
update_recipe_mode_for_pedigree(
    {
        active_parameter_href => \%active_parameter,
        recipes_ref           => [qw{ blobfish }],
        sample_info_href      => \%sample_info,
    }
);

## Then turn off blobfish
is( $active_parameter{blobfish}, 0, q{Recipe is not compatible with pedigree} );

## Given an active recipe
$active_parameter{blobfish} = 1;
delete $sample_info{sample}{ADM1059A1};
delete $sample_info{sample}{ADM1059A2};

## When single sample in pedigree
update_recipe_mode_for_pedigree(
    {
        active_parameter_href => \%active_parameter,
        recipes_ref           => [qw{ blobfish }],
        sample_info_href      => \%sample_info,
    }
);

## Then turn off blobfish
is( $active_parameter{blobfish}, 0, q{Recipe is not compatible with pedigree} );
done_testing();
