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
use MIP::Constants qw{ $COMMA $EMPTY_STR $SPACE $UNDERSCORE };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Processmanagement::Processes} => [qw{ add_to_job_id_dependency_string }],

    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Processmanagement::Processes qw{ add_to_job_id_dependency_string };

diag(   q{Test add_to_job_id_dependency_string from Processes.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Base arguments
my $sample_id         = q{sample2};
my $path              = q{MAIN};
my $case_id_chain_key = q{case1} . $UNDERSCORE . $path;

my %job_id = (
    $case_id_chain_key => {
        q{sample1} . $UNDERSCORE . $path => [qw{job_id_1 job_id_2}],
        q{sample2} . $UNDERSCORE . $path => [qw{job_id_3}],
        q{sample3} . $UNDERSCORE . $path => [qw{job_id_4 job_id_5 job_id_8}],
        q{sample4} . $UNDERSCORE . $path => [undef],
        $case_id_chain_key               => [qw{job_id_6}],
    },
);

### Sample job

## Add 1 job_id to job_id_string
my $sample_id_chain_key = $sample_id . $UNDERSCORE . $path;

my $job_ids_string = add_to_job_id_dependency_string(
    {
        case_id_chain_key => $case_id_chain_key,
        chain_key         => $sample_id_chain_key,
        job_id_href       => \%job_id,
    }
);
my $expected_job_id_string = q{:job_id_3};
is( $job_ids_string, $expected_job_id_string, q{Added 1 job_id to job_id_string} );

## Add 2 job_ids to job_id_string
$sample_id           = q{sample1};
$sample_id_chain_key = $sample_id . $UNDERSCORE . $path;

## Add to job_id string
$job_ids_string = add_to_job_id_dependency_string(
    {
        case_id_chain_key => $case_id_chain_key,
        chain_key         => $sample_id_chain_key,
        job_id_href       => \%job_id,
    }
);

$expected_job_id_string = q{:job_id_1:job_id_2};

is( $job_ids_string, $expected_job_id_string, q{Added 2 job_ids to job_id_string} );

## Add 3 job_ids to job_id_string
$sample_id           = q{sample3};
$sample_id_chain_key = $sample_id . $UNDERSCORE . $path;

## Add to job_id string
$job_ids_string = add_to_job_id_dependency_string(
    {
        job_id_href       => \%job_id,
        case_id_chain_key => $case_id_chain_key,
        chain_key         => $sample_id_chain_key,
    }
);

$expected_job_id_string = q{:job_id_4:job_id_5:job_id_8};

is( $job_ids_string, $expected_job_id_string, q{Added 3 job_ids to job_id_string} );

## Do not add undef job_ids to job_id_string
$sample_id           = q{sample4};
$sample_id_chain_key = $sample_id . $UNDERSCORE . $path;

## Add to job_id string
$job_ids_string = add_to_job_id_dependency_string(
    {
        job_id_href       => \%job_id,
        case_id_chain_key => $case_id_chain_key,
        chain_key         => $sample_id_chain_key,
    }
);

$expected_job_id_string = $EMPTY_STR;

is( $job_ids_string, $expected_job_id_string, q{Nothing was added to job_id_string} );

done_testing();
