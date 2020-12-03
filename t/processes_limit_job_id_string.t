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
use Readonly;

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
        q{MIP::Processmanagement::Processes} => [qw{ limit_job_id_string }],
        q{MIP::Test::Fixtures}               => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Processmanagement::Processes qw{ limit_job_id_string };

diag(   q{Test limit_job_id_string from Processes.pm v}
      . $MIP::Processmanagement::Processes::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

use MIP::Processmanagement::Processes qw{ limit_job_id_string };

## Constants
Readonly my $MAX_JOB_IDS_TO_TRACK      => q{1001};
Readonly my $OVER_MAX_JOB_IDS_TO_TRACK => q{1200};

# Create job_ids array
my @job_ids = ( 0 .. $OVER_MAX_JOB_IDS_TO_TRACK );

## Base arguments
my $case_id             = q{case1};
my $sample_id           = q{sample1};
my $path                = q{MAIN};
my $case_id_chain_key   = $case_id . $UNDERSCORE . $path;
my $sample_id_chain_key = $sample_id . $UNDERSCORE . $path;
my $pan_chain_key       = $case_id_chain_key . $UNDERSCORE . $sample_id_chain_key;

my %job_id = (
    $case_id_chain_key => {
        q{sample1} . $UNDERSCORE . $path => [qw{job_id_1 job_id_2}],
        q{sample2} . $UNDERSCORE . $path => [qw{job_id_3}],
        q{sample3} . $UNDERSCORE . $path => [qw{job_id_4 job_id_5 job_id_8}],
        q{sample4} . $UNDERSCORE . $path => [undef],
        $pan_chain_key                   => [qw{job_id_1 job_id_2}],
        $case_id_chain_key               => [qw{job_id_6}],
    },
    q{ALL} => { q{ALL} => [@job_ids], }
);

## When reducing the size of the job_ids array
limit_job_id_string(
    {
        job_id_href => \%job_id,
    }
);

my $result_ref      = scalar @{ $job_id{q{ALL}}{q{ALL}} };
my $expected_result = $MAX_JOB_IDS_TO_TRACK;

## Then the number of job ids is reduced
is( $result_ref, $expected_result, q{Limited nr of job_ids in job_id chain} );

## When adding job_ids from case MAIN chain to job_id_string
limit_job_id_string(
    {
        job_id_href       => \%job_id,
        case_id_chain_key => $case_id_chain_key,
        chain_key         => $sample_id_chain_key,
    }
);

$result_ref      = scalar @{ $job_id{$case_id_chain_key}{$sample_id_chain_key} };
$expected_result = q{2};

## Then two jobs should be returned
is( $result_ref, $expected_result, q{Keept job_ids in job_id chain} );

done_testing();
