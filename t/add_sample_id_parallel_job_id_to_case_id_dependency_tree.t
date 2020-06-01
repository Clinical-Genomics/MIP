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
use MIP::Constants qw{ $COMMA $SPACE $UNDERSCORE };
use MIP::Test::Fixtures qw{ test_standard_cli };

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
        q{MIP::Processmanagement::Processes} => [qw{ add_sample_id_parallel_job_id_to_case_id_dependency_tree }],
        q{MIP::Test::Fixtures}   => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Processmanagement::Processes qw{ add_sample_id_parallel_job_id_to_case_id_dependency_tree };

diag(   q{Test add_sample_id_parallel_job_id_to_case_id_dependency_tree from Processes.pm v}
      . $MIP::Processmanagement::Processes::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my $case_id             = q{case1};
my $path                = q{MAIN};
my $sample_id           = q{sample2};

my $case_id_chain_key   = $case_id . $UNDERSCORE . $path;
my $sample_id_chain_key = $sample_id . $UNDERSCORE . $path;
my $parallel_processes_index        = 0;
my $sample_id_parallel_chain_key =
  $sample_id . $UNDERSCORE . q{parallel} . $UNDERSCORE . $path . $parallel_processes_index;

my %max_parallel_processes_count = (
    sample1 => 1,
    sample2 => 0,
    sample3 => 1,
);
my %job_id = (
    $case_id_chain_key => {
        $sample_id_chain_key          => [qw{job_id_1 job_id_2}],
        q{sample2_MAIN}               => [qw{job_id_3}],
        q{sample3_MAIN}               => [qw{job_id_4 job_id_5}],
        $sample_id_parallel_chain_key => [qw{job_id_10 job_id_11}],
        $case_id_chain_key            => [qw{job_id_6}],
    },
);

## Given a sample_id
## When no parallel jobs

add_sample_id_parallel_job_id_to_case_id_dependency_tree(
    {
        case_id_chain_key       => $case_id_chain_key,
        job_id_href             => \%job_id,
        max_parallel_processes_count_href => \%max_parallel_processes_count,
        path                    => $path,
        sample_id               => $sample_id . q{no_parallel},
    }
);

my $case_id_job_str = join $SPACE, @{ $job_id{$case_id_chain_key}{$case_id_chain_key} };

## Then return only case id job_ids and no sample id jobs
is( $case_id_job_str, q{job_id_6}, q{No sample_id parallel job_ids} );

## When a sample has parallel jobs
add_sample_id_parallel_job_id_to_case_id_dependency_tree(
    {
        case_id_chain_key       => $case_id_chain_key,
        job_id_href             => \%job_id,
        max_parallel_processes_count_href => \%max_parallel_processes_count,
        path                    => $path,
        sample_id               => $sample_id,
    }
);

my $sample_id_parallel_push_result = join $SPACE,
  @{ $job_id{$case_id_chain_key}{$case_id_chain_key} };

## Then return case id job_ids and sample id parallel jobs
is(
    $sample_id_parallel_push_result,
    q{job_id_6 job_id_10 job_id_11},
    q{Pushed parallel job_id to sample_id}
);

done_testing();
