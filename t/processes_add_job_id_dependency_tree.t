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
        q{MIP::Processmanagement::Processes} => [qw{ add_job_id_dependency_tree }],

    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Processmanagement::Processes qw{ add_job_id_dependency_tree };

diag(   q{Test add_job_id_dependency_tree from Processes.pm}
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

my %job_id = (
    $case_id_chain_key => {
        $sample_id_chain_key => [qw{job_id_1 job_id_2}],
        q{sample2_MAIN}      => [qw{job_id_3}],
        q{sample3_MAIN}      => [qw{job_id_4 job_id_5}],
        $case_id_chain_key   => [qw{job_id_6}],
    },
);

### Add to arbitrary chain key

my $job_id_returned = q{job_id_7};

add_job_id_dependency_tree(
    {
        chain_key       => $case_id_chain_key,
        job_id_href     => \%job_id,
        job_id_returned => $job_id_returned,
    }
);

my $result = join $SPACE, @{ $job_id{$case_id_chain_key}{$case_id_chain_key} };
is( $result, q{job_id_6 job_id_7}, q{Pushed to job_id to arbitrary chain key} );

done_testing();
