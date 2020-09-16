#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $COLON $EMPTY_STR $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

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
        q{MIP::Processmanagement::Processes} => [qw{ track_job_id_status }],
        q{MIP::Test::Fixtures}               => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Processmanagement::Processes qw{ track_job_id_status };

diag(   q{Test track_job_id_status from Processes.pm v}
      . $MIP::Processmanagement::Processes::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a filehandle
# For storing info to write
my $file_content;

## Store file content in memory by using referenced variable
open my $filehandle, q{>}, \$file_content
  or croak q{Cannot write to} . $SPACE . $file_content . $COLON . $SPACE . $OS_ERROR;

## When no log_file_path and job_ids
my $return = track_job_id_status(
    {
        filehandle => $filehandle,
    }
);
is( $return, 0, q{Do not write command to track job_ids} );

## Given a job_id and a log path
my @job_ids       = (qw{ job_id_test });
my $log_file_path = catfile( $Bin, qw{ data test_data log_test} );

## When no submission profile
$return = track_job_id_status(
    {
        filehandle    => $filehandle,
        job_ids_ref   => \@job_ids,
        log_file_path => $log_file_path,
    }
);
is( $return, undef, q{Do not write command to track job_ids} );

## Given a job_id and a log path

## When using a slurm submission profile
my @sacct_format_fields = qw{
  jobid jobname%50 account partition alloccpus TotalCPU elapsed start end state exitcode };

track_job_id_status(
    {
        filehandle              => $filehandle,
        job_ids_ref             => \@job_ids,
        log_file_path           => $log_file_path,
        submission_profile      => q{slurm},
        sacct_format_fields_ref => \@sacct_format_fields,
    }
);

## Close the filehandle
close $filehandle;

## Then track progress cmd should be written to file
my ($wrote_string) = $file_content =~ /(sacct \s+ --format=)/xms;

ok( $wrote_string, q{Wrote instruction to track slurm job_ids progress} );

done_testing();
