#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Path qw{ remove_tree };
use File::Spec::Functions qw{ catdir catfile };
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
use MIP::Constants qw{ $COMMA $DOT $SPACE $UNDERSCORE };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

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
        q{MIP::Processmanagement::Processes} => [qw{ write_job_ids_to_file }],
        q{MIP::Test::Fixtures}               => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Processmanagement::Processes qw{ write_job_ids_to_file };

diag(   q{Test write_job_ids_to_file from Processes.pm v}
      . $MIP::Processmanagement::Processes::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { no_screen => 1, } );

## Given a log file when no job ids
my $case_id          = q{case_1};
my %active_parameter = (
    case_id  => $case_id,
    log_file => catfile( $Bin, q{a_test.log} ),
);
my $date_time_stamp = q{1999-12-31T12:00:00};
my %job_id;
my $is_ok = write_job_ids_to_file(
    {
        active_parameter_href => \%active_parameter,
        date_time_stamp       => $date_time_stamp,
        job_id_href           => \%job_id,
    }
);

## Then return false
is( $is_ok, undef, q{Skip when no job ids} );

## Given a log file when job_ids
%job_id = ( ALL => { ALL => [ qw{ job_id_1 job_id_2 }, undef, ], }, );
$is_ok  = write_job_ids_to_file(
    {
        active_parameter_href => \%active_parameter,
        date_time_stamp       => $date_time_stamp,
        job_id_href           => \%job_id,
    }
);

## Then return true
ok( $is_ok, q{Wrote job ids file for job ids} );

## Then job_id file should exist
my $log_dir = dirname( $active_parameter{log_file} );
my $job_ids_file =
  catfile( $log_dir, q{slurm_job_ids} . $UNDERSCORE . $date_time_stamp . $DOT . q{yaml} );
ok( -e $job_ids_file, q{Created file} );
## Clean-up
remove_tree($job_ids_file);

done_testing();
