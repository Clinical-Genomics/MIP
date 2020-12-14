#! /usr/bin/env perl

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
use MIP::Constants qw{ $COMMA $SPACE $UNDERSCORE };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Processmanagement::Processes} =>
          [qw{ add_parallel_chains_job_ids_to_job_id_dependency_string }],

    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Processmanagement::Processes qw{ add_parallel_chains_job_ids_to_job_id_dependency_string };

diag(   q{Test add_parallel_chains_job_ids_to_job_id_dependency_string from Processes.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my $case_id             = q{case1};
my $sample_id           = q{sample2};
my $path                = q{MAIN};
my $case_id_chain_key   = $case_id . $UNDERSCORE . $path;
my $sample_id_chain_key = $sample_id . $UNDERSCORE . $path;
my $infile_index        = 0;
my $sample_id_parallel_chain_key =
  $sample_id . $UNDERSCORE . q{parallel} . $UNDERSCORE . $path . $infile_index;

my $case_id_chain_key_other = $case_id . $UNDERSCORE . q{OTHER};

my @sample_ids      = ( $sample_id, qw{sample4} );
my @parallel_chains = qw{OTHER};

my %infile_lane_prefix = (
    sample1 => [qw{1_lane1 1_lane2}],
    sample2 => [qw{2_lane1}],
    sample3 => [qw{3_lane4 3_lane5}],
);

my %job_id = (
    $case_id_chain_key => {
        $sample_id_chain_key          => [qw{job_id_1 job_id_2}],
        q{sample2_MAIN}               => [qw{job_id_3}],
        q{sample3_MAIN}               => [qw{job_id_4 job_id_5}],
        $sample_id_parallel_chain_key => [qw{job_id_10 job_id_11}],
        $case_id_chain_key            => [qw{job_id_6}],
    },
    $case_id_chain_key_other => {
        q{sample4_OTHER}         => [qw{job_id_12 job_id_13}],
        $case_id_chain_key_other => [qw{job_id_14}],
    },
);

### Parallel chain jobs

my $job_id_string = add_parallel_chains_job_ids_to_job_id_dependency_string(
    {
        case_id             => $case_id,
        job_id_href         => \%job_id,
        parallel_chains_ref => \@parallel_chains,
        sample_ids_ref      => \@sample_ids,
    }
);

is(
    $job_id_string,
    q{:job_id_12:job_id_13:job_id_14},
    q{Added parallel chains job_ids to job_id string}
);

done_testing();
