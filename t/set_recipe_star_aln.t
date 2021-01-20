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
use MIP::Test::Fixtures qw{ test_mip_hashes };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Analysis}       => [qw{ set_recipe_star_aln }],
        q{MIP::Test::Fixtures} => [qw{ test_mip_hashes }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Analysis::Star_aln qw{ analysis_star_aln analysis_star_aln_mixed };
use MIP::Analysis qw{ set_recipe_star_aln };

diag(   q{Test set_recipe_star_aln from Analysis.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given an analysis recipe and a mip_file_format
my %analysis_recipe;
my $mip_file_format = q{ADM1059A1_161011_HHJJCCCXY_NAATGCGC_lane1};

## Given a sample with file info
my %file_info = test_mip_hashes( { mip_hash_name => q{file_info}, } );

## Given a sample id
my $sample_id  = q{ADM1059A1};
my @sample_ids = ( $sample_id, );

## When setting star_aln recipe
set_recipe_star_aln(
    {
        analysis_recipe_href => \%analysis_recipe,
        file_info_href       => \%file_info,
        sample_ids_ref       => \@sample_ids,
    }
);

## Then use analysis_star_aln recipe
is( $analysis_recipe{star_aln}, \&analysis_star_aln, q{Set star_aln single sequence run type} );

## Given a sample with mixed single and paired end fastq files
push @{ $file_info{$sample_id}{no_direction_infile_prefixes} }, $mip_file_format;
$file_info{$sample_id}{$mip_file_format}{sequence_run_type} = q{paired-end};

## When setting star_aln recipe
set_recipe_star_aln(
    {
        analysis_recipe_href => \%analysis_recipe,
        file_info_href       => \%file_info,
        sample_ids_ref       => \@sample_ids,
    }
);

## Then use analysis_star_aln mixed recipe
is( $analysis_recipe{star_aln},
    \&analysis_star_aln_mixed, q{Set star_aln multiple sequence run type} );

done_testing();
