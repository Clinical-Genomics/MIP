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
use MIP::Test::Fixtures qw{ test_mip_hashes };

## Constants
Readonly my $COMMA      => q{,};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Processmanagement::Processes} =>
          [qw{ create_job_id_string_for_sample_id }],
        q{MIP::Test::Fixtures} => [qw{ test_mip_hashes }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Processmanagement::Processes qw{ create_job_id_string_for_sample_id };

diag(   q{Test create_job_id_string_for_sample_id from Processes.pm}
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
my $sample_id_chain_key   = $sample_id . $UNDERSCORE . $path;
my $sample_id_parallel_chain_key =
  $sample_id . $UNDERSCORE . q{parallel} . $UNDERSCORE . $path . $sbatch_script_tracker;

my %job_id = test_mip_hashes( { mip_hash_name => q{job_id}, } );

## Given job ids from MAIN chain, when using sample id

## Add job_ids from MAIN chain to job_id_string
my $job_ids_string = create_job_id_string_for_sample_id(
    {
        case_id             => $case_id,
        case_id_chain_key   => $case_id_chain_key,
        job_id_href         => \%job_id,
        path                => $path,
        sample_id           => $sample_id,
        sample_id_chain_key => $sample_id_chain_key,
    }
);

## Then add job_ids for sample1 from MAIN
my $expected_job_id_string = q{:job_id_1:job_id_2};
is( $job_ids_string, $expected_job_id_string, q{Added job_id from MAIN job_id chain} );

## Given job id string using job ids from other chain with no previous job ids for sample id
my $path_other                = q{other};
my $case_id_chain_key_other   = $case_id . $UNDERSCORE . $path_other;
my $sample_id_chain_key_other = $sample_id . $UNDERSCORE . $path_other;

## Add job_ids from MAIN chain to job_id_string
$job_ids_string = create_job_id_string_for_sample_id(
    {
        case_id             => $case_id,
        case_id_chain_key   => $case_id_chain_key_other,
        job_id_href         => \%job_id,
        path                => $path_other,
        sample_id           => $sample_id,
        sample_id_chain_key => $sample_id_chain_key_other,
    }
);

## Then add job_ids for sample1 inherited from MAIN chain job ids
my $expected_job_id_string_empty_other = q{:job_id_1:job_id_2};
is(
    $job_ids_string,
    $expected_job_id_string_empty_other,
    q{Added job_id from MAIN job_id chain initializing other chain}
);

### Inherit job_ids from other chain

%job_id = (
    $case_id_chain_key_other => {
        q{sample1} . $UNDERSCORE . $path_other => [qw{job_id_9 job_id_10}],
    },
);

## Given job id string using job ids from other chain for sample id
$job_ids_string = create_job_id_string_for_sample_id(
    {
        case_id             => $case_id,
        case_id_chain_key   => $case_id_chain_key_other,
        job_id_href         => \%job_id,
        path                => $path_other,
        sample_id           => $sample_id,
        sample_id_chain_key => $sample_id_chain_key_other,
    }
);

my $expected_job_id_string_other = q{:job_id_9:job_id_10};

## Then add job_ids for other chain
is( $job_ids_string, $expected_job_id_string_other,
    q{Added job_id from other job_id chain} );

## Given job id string using parallel job ids from other chain with no previous job ids for sample id

## Clean-up for new test
%job_id = ();

$job_id{$case_id_chain_key}{$sample_id_parallel_chain_key} =
  [qw{ job_id_11 job_id_12 }];

## Add job_ids from MAIN chain to job_id_string
$job_ids_string = create_job_id_string_for_sample_id(
    {
        case_id               => $case_id,
        case_id_chain_key     => $case_id_chain_key_other,
        job_id_href           => \%job_id,
        path                  => $path_other,
        sample_id             => $sample_id,
        sample_id_chain_key   => $sample_id_chain_key_other,
        sbatch_script_tracker => $sbatch_script_tracker,
    }
);

## Then add job_ids for sample1 inherited from MAIN chain parallel job ids
my $expected_job_id_string_parallel_other = q{:job_id_11:job_id_12};
is(
    $job_ids_string,
    $expected_job_id_string_parallel_other,
    q{Added parallel job_id from MAIN job_id chain initializing other chain}
);

done_testing();
