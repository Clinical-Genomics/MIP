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
          [qw{ add_parallel_job_ids_to_job_id_dependency_string }],

    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Processmanagement::Processes qw{ add_parallel_job_ids_to_job_id_dependency_string };

diag(   q{Test add_parallel_job_ids_to_job_id_dependency_string from Processes.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my $case_id               = q{case1};
my $sample_id             = q{sample2};
my $path                  = q{MAIN};
my $case_id_chain_key     = $case_id . $UNDERSCORE . $path;
my $sample_id_chain_key   = $sample_id . $UNDERSCORE . $path;
my $sbatch_script_tracker = 0;
my $case_id_parallel_chain_key =
  $case_id . $UNDERSCORE . q{parallel} . $UNDERSCORE . $path . $sbatch_script_tracker;

my %job_id = (
    $case_id_chain_key => {
        $sample_id_chain_key        => [qw{job_id_1 job_id_2}],
        q{sample2_MAIN}             => [qw{job_id_3}],
        q{sample3_MAIN}             => [qw{job_id_4 job_id_5}],
        $case_id_parallel_chain_key => [qw{job_id_10 job_id_11}],
        $case_id_chain_key          => [qw{job_id_6}],
    },
);

### Parallel jobs

my $no_job_id_string = add_parallel_job_ids_to_job_id_dependency_string(
    {
        case_id_chain_key => $case_id_chain_key . q{no_parallel},
        job_id_href       => \%job_id,
    }
);

is( $no_job_id_string, undef, q{No parallel job_ids} );

my $job_id_string = add_parallel_job_ids_to_job_id_dependency_string(
    {
        case_id_chain_key => $case_id_chain_key,
        job_id_href       => \%job_id,
    }
);

is( $job_id_string, q{:job_id_10:job_id_11}, q{Added parallel job_ids to job id string} );

done_testing();
