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
our $VERSION = 1.02;

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
        q{MIP::Processmanagement::Processes} =>
          [qw{ add_parallel_job_id_to_sample_id_dependency_tree }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Processmanagement::Processes
  qw{add_parallel_job_id_to_sample_id_dependency_tree};

diag(   q{Test add_parallel_job_id_to_sample_id_dependency_tree from Processes.pm v}
      . $MIP::Processmanagement::Processes::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my $sample_id           = q{sample1};
my $path                = q{MAIN};
my $case_id_chain_key   = q{case1} . $UNDERSCORE . $path;
my $sample_id_chain_key = $sample_id . $UNDERSCORE . $path;
my $chain_key_parallel  = q{parallel};

my %max_parallel_processes_count = (
    sample1 => 1,
    sample2 => 0,
    sample3 => 1,
);
my %job_id = (
    $case_id_chain_key => {
        $sample_id_chain_key => [qw{job_id_1 job_id_2}],
        q{sample2_MAIN}      => [qw{job_id_3}],
        q{sample3_MAIN}      => [qw{job_id_4 job_id_5}],
        $case_id_chain_key   => [qw{job_id_6}],
    },
);

### No previous parallel job_ids
add_parallel_job_id_to_sample_id_dependency_tree(
    {
        case_id_chain_key                 => $case_id_chain_key,
        job_id_href                       => \%job_id,
        max_parallel_processes_count_href => \%max_parallel_processes_count,
        path                              => $path,
        sample_id                         => $sample_id,
        sample_id_chain_key               => $sample_id_chain_key,
    }
);
my $no_parallel_push_result = join $SPACE,
  @{ $job_id{$case_id_chain_key}{$sample_id_chain_key} };
is( $no_parallel_push_result, q{job_id_1 job_id_2}, q{No parallel job_id} );

### Previous parallel jobs
## Set-up previous parallel job
PARALLEL_PROCESS:
foreach my $parallel_processes_index ( 0 .. $max_parallel_processes_count{$sample_id} ) {

    # Set key
    my $chain_key_parallel_job =
        $sample_id
      . $UNDERSCORE
      . $chain_key_parallel
      . $UNDERSCORE
      . $path
      . $parallel_processes_index;

    push @{ $job_id{$case_id_chain_key}{$chain_key_parallel_job} },
      q{job_id_} . $parallel_processes_index;
}

add_parallel_job_id_to_sample_id_dependency_tree(
    {
        case_id_chain_key                 => $case_id_chain_key,
        job_id_href                       => \%job_id,
        max_parallel_processes_count_href => \%max_parallel_processes_count,
        path                              => $path,
        sample_id                         => $sample_id,
        sample_id_chain_key               => $sample_id_chain_key,
    }
);

my $parallel_push_result = join $SPACE,
  @{ $job_id{$case_id_chain_key}{$sample_id_chain_key} };
is(
    $parallel_push_result,
    q{job_id_1 job_id_2 job_id_0 job_id_1},
    q{Push parallel job_id}
);

done_testing();
