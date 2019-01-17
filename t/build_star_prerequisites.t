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
use Modern::Perl qw{ 2014 };
use Readonly;
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_log test_mip_hashes test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $COLON => q{:};
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Recipes::Build::Star_prerequisites} => [qw{ build_star_prerequisites }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Build::Star_prerequisites qw{ build_star_prerequisites };

diag(   q{Test build_star_prerequisites from Star_prerequisites.pm v}
      . $MIP::Recipes::Build::Star_prerequisites::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log();

## Given build parameters
my $recipe_name = q{star_fusion};

my %active_parameter = test_mip_hashes(
    {
        mip_hash_name => q{active_parameter},
        recipe_name   => $recipe_name,
    }
);
$active_parameter{star_fusion}               = 1;
$active_parameter{star_aln_reference_genome} = q{_star_genome_dir};
my %file_info = test_mip_hashes(
    {
        mip_hash_name => q{file_info},
        recipe_name   => $recipe_name,
    }
);
my %infile_lane_prefix;
my %job_id;
my %parameter = test_mip_hashes( { mip_hash_name => q{recipe_parameter}, } );

my %sample_info;

# Special case
$active_parameter{transcript_annotation} = q{GRCH37_transcripts.gtf};

trap {
    build_star_prerequisites(
        {
            active_parameter_href        => \%active_parameter,
            file_info_href               => \%file_info,
            infile_lane_prefix_href      => \%infile_lane_prefix,
            job_id_href                  => \%job_id,
            log                          => $log,
            parameter_href               => \%parameter,
            parameter_build_suffixes_ref => \@{ $file_info{star_aln_reference_genome} },
            recipe_name                  => $recipe_name,
            sample_info_href             => \%sample_info,
        }
    )
};

## Then broadcast info log message
my $log_msg = q{Will\s+try\s+to\s+create\s+required\s+human_genome.fasta\s+star\s+files};
like( $trap->stderr, qr/$log_msg/msx, q{Broadcast star_fusion log message} );

done_testing();
