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
        q{MIP::Processmanagement::Processes} => [qw{ add_sample_ids_job_ids_to_job_id_dependency_string }],
        q{MIP::Test::Fixtures}   => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Processmanagement::Processes qw{ add_sample_ids_job_ids_to_job_id_dependency_string };

diag(   q{Test add_sample_ids_job_ids_to_job_id_dependency_string from Processes.pm v}
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
my @samples;
my $sample_id           = q{sample2};
my $parallel_processes_index        = 0;

my $case_id_chain_key   = q{case1} . $UNDERSCORE . $path;
my $sample_id_parallel_chain_key =
  $sample_id . $UNDERSCORE . q{parallel} . $UNDERSCORE . $path . $parallel_processes_index;

my %max_parallel_processes_count = (
    sample1 => 1,
    sample2 => 0,
    sample3 => 1,
);
my %job_id = (
    $case_id_chain_key => {
        q{sample1} . $UNDERSCORE . $path => [qw{job_id_1 job_id_2}],
        q{sample2} . $UNDERSCORE . $path => [qw{job_id_3}],
        q{sample3} . $UNDERSCORE . $path => [qw{job_id_4 job_id_5 job_id_8}],
        q{sample4} . $UNDERSCORE . $path => [undef],
        $sample_id_parallel_chain_key    => [qw{job_id_10 job_id_11}],
        $case_id_chain_key               => [qw{job_id_6}],
    },
);

## Given a no samples

## When sample has no parallel jobs

my $job_ids_string = add_sample_ids_job_ids_to_job_id_dependency_string(
    {
        case_id                 => $case_id,
        case_id_chain_key       => $case_id_chain_key,
        job_id_href             => \%job_id,
        max_parallel_processes_count_href => \%max_parallel_processes_count,
        path                    => $path,
        sample_ids_ref          => \@samples,
    }
);
my $expected_job_id_string;

## Then return undef job_ids_string
is( $job_ids_string, $expected_job_id_string,
    q{No sample_id and no parallel jobs to job_id_string} );

## Given one sample
@samples = qw{sample1};

##  When no sample_id previous parallel jobs exists, but two parallel processes are started
$job_ids_string = add_sample_ids_job_ids_to_job_id_dependency_string(
    {
        case_id                 => $case_id,
        case_id_chain_key       => $case_id_chain_key,
        job_id_href             => \%job_id,
        max_parallel_processes_count_href => \%max_parallel_processes_count,
        path                    => $path,
        sample_ids_ref          => \@samples,
    }
);
$expected_job_id_string = q{:job_id_1:job_id_2};

## Then return job_ids_string for sample with new parallel processes
is( $job_ids_string, $expected_job_id_string,
    q{Added sample_id and no parallel jobs to job_id_string} );

## Given three samples
@samples = qw{sample1 sample2 sample3};

## When previous parallel jobs exist
$job_ids_string = add_sample_ids_job_ids_to_job_id_dependency_string(
    {
        case_id                 => $case_id,
        case_id_chain_key       => $case_id_chain_key,
        job_id_href             => \%job_id,
        max_parallel_processes_count_href => \%max_parallel_processes_count,
        path                    => $path,
        sample_ids_ref          => \@samples,
    }
);
$expected_job_id_string =
  q{:job_id_1:job_id_2:job_id_3:job_id_10:job_id_11:job_id_4:job_id_5:job_id_8};

## Then return job_ids for all samples as well as existing previous parallel job_ids and newly started per sample
is( $job_ids_string, $expected_job_id_string,
    q{Added sample_ids and parallel jobs to job_id_string} );

done_testing();

