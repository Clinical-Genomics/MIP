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
        q{MIP::Set::Analysis}  => [qw{ set_recipe_star_aln }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Analysis::Star_aln qw{ analysis_star_aln analysis_star_aln_mixed };
use MIP::Set::Analysis qw{ set_recipe_star_aln };

diag(   q{Test set_recipe_star_aln from Analysis.pm v}
      . $MIP::Set::Analysis::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a sample with only paired end fastq files.
my $sample_id               = q{test};
my %infile_lane_prefix_hash = ( $sample_id => [qw{ lane1 lane2 }], );
my %sample_info             = (
    sample => {
        $sample_id => {
            file => {
                lane1 => {
                    sequence_run_type => q{paired-end},
                },
                lane2 => {
                    sequence_run_type => q{paired-end},
                },
            },
        },
    },
);
my %analysis_recipe;

## Then use analysis_star_aln recipe
set_recipe_star_aln(
    {
        analysis_recipe_href    => \%analysis_recipe,
        infile_lane_prefix_href => \%infile_lane_prefix_hash,
        sample_info_href        => \%sample_info,
    }
);
is( $analysis_recipe{star_aln},
    \&analysis_star_aln, q{Set star_aln single sequence run type} );

## Given a sample with mixed single and paired end fastq files.
$sample_info{sample}{$sample_id}{file}{lane2}{sequence_run_type} = q{single-end};

## Then use analysis_star_aln recipe_mixed
set_recipe_star_aln(
    {
        analysis_recipe_href    => \%analysis_recipe,
        infile_lane_prefix_href => \%infile_lane_prefix_hash,
        sample_info_href        => \%sample_info,
    }
);
is( $analysis_recipe{star_aln},
    \&analysis_star_aln_mixed, q{Set star_aln multiple sequence run type} );

done_testing();
