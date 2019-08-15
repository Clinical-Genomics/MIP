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
use MIP::Test::Fixtures qw{ test_mip_hashes test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $COMMA      => q{,};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Processmanagement::Processes} => [qw{ create_job_id_string_for_case_id }],
        q{MIP::Test::Fixtures}               => [qw{ test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Processmanagement::Processes qw{ create_job_id_string_for_case_id };

diag(   q{Test create_job_id_string_for_case_id from Processes.pm v}
      . $MIP::Processmanagement::Processes::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my $case_id               = q{case1};
my $path                  = q{MAIN};
my $sample_id             = q{sample1};
my $sbatch_script_tracker = 0;
my $case_id_chain_key     = $case_id . $UNDERSCORE . $path;
my $case_id_parallel_chain_key =
  $case_id . $UNDERSCORE . q{parallel} . $UNDERSCORE . $path . $sbatch_script_tracker;
my $parallel_chain_key = q{parallel} . $UNDERSCORE . $path . $sbatch_script_tracker;
my @parallel_chains    = ($parallel_chain_key);
my @sample_ids         = ($sample_id);
my %infile_lane_prefix;
my %job_id = test_mip_hashes( { mip_hash_name => q{job_id}, } );

## Given job ids from MAIN chain, when existing job id for case id
my $job_ids_string = create_job_id_string_for_case_id(
    {
        case_id                 => $case_id,
        case_id_chain_key       => $case_id_chain_key,
        infile_lane_prefix_href => \%infile_lane_prefix,
        job_id_href             => \%job_id,
        path                    => $path,
        sample_ids_ref          => \@sample_ids,
    }
);

## Then add job_ids for case1 from MAIN
my $expected_job_id_string = q{:job_id_6};
is( $job_ids_string, $expected_job_id_string, q{Added job_id from MAIN job_id chain} );

## Given job id string using job ids from other chain with no previous job ids for case id
my $path_other              = q{other};
my $case_id_chain_key_other = $case_id . $UNDERSCORE . $path_other;
my $job_ids_string_other    = create_job_id_string_for_case_id(
    {
        case_id                 => $case_id,
        case_id_chain_key       => $case_id_chain_key_other,
        infile_lane_prefix_href => \%infile_lane_prefix,
        job_id_href             => \%job_id,
        path                    => $path_other,
        sample_ids_ref          => \@sample_ids,
    }
);

## Then add job_ids for case1 from MAIN
my $expected_job_id_string_other = q{:job_id_6};
is( $job_ids_string_other, $expected_job_id_string_other,
    q{Added job_id from to other job_id chain} );

## Given job id string using parallel job ids from other chain with a previous job ids for case id
$job_id{$case_id_parallel_chain_key}{$case_id_parallel_chain_key} =
  [qw{ job_id_11 job_id_12 }];

my $job_ids_string_other_parallel = create_job_id_string_for_case_id(
    {
        case_id                 => $case_id,
        case_id_chain_key       => $case_id_chain_key_other,
        infile_lane_prefix_href => \%infile_lane_prefix,
        job_id_href             => \%job_id,
        parallel_chains_ref     => \@parallel_chains,
        path                    => $path_other,
        sample_ids_ref          => \@sample_ids,
    }
);

## Then add job_ids for case1 from MAIN
my $expected_job_id_string_other_parallel = q{:job_id_6:job_id_11:job_id_12};
is(
    $job_ids_string_other_parallel,
    $expected_job_id_string_other_parallel,
    q{Added job_id from main and main parallel to other job_id chain}
);

done_testing();
