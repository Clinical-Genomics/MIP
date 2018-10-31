#!/usr/bin/env perl

use 5.026;
use Carp;
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use IPC::Cmd qw(can_run run);
use Params::Check qw{ allow check last_error };
use Test::More;
use charnames qw{ :full :short };
use open qw{ :encoding(UTF-8) :std };
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2014 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.01;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $COMMA   => q{,};
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };
Readonly my $TAB     => qq{\t};

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Workloadmanager::Slurm} => [qw{ slurm_reformat_sacct_output}],
        q{MIP::Test::Fixtures}         => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Language::Shell qw{ build_shebang };
use MIP::Gnu::Bash qw{ gnu_set };
use MIP::Workloadmanager::Slurm qw{ slurm_reformat_sacct_output };

diag(   q{Test slurm_reformat_sacct_output from Slurm.pm v}
      . $MIP::Workloadmanager::Slurm::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

###  Set up test data parameters

## Test input from sacct command
my $test_data_file_path =
  catfile( $Bin, qw{ data test_data slurm_reformat_sacct_output_input.txt } );

## Header for output status file
my @reformat_sacct_headers =
  qw{ JobID JobName Account Partition AllocCPUS TotalCPU Elapsed Start End State ExitCode };

## Write test bash script to this file
my $bash_file_path = catfile( $Bin, q{test_slurm_reformat_sacct_output.sh} );

## File wit sacct output to reformat
my $log_file_path = catfile( $Bin, q{slurm_reformat_sacct_output.txt} );

## Output file with reformated sacct output
my $outfile_path = $log_file_path . q{.status};

### Start test

## Create anonymous filehandle
my $FILEHANDLE = IO::Handle->new();

## Open filehandle
open $FILEHANDLE, q{>}, $bash_file_path
  or croak( q{Cannot write to '} . $bash_file_path . q{' :} . $OS_ERROR . $NEWLINE );

## Write reformat  commando to bash file
my @commands = ( q{less}, $test_data_file_path );
_build_test_file_recipe(
    {
        bash_file_path             => $bash_file_path,
        commands_ref               => \@commands,
        FILEHANDLE                 => $FILEHANDLE,
        log_file_path              => $log_file_path,
        reformat_sacct_headers_ref => \@reformat_sacct_headers,
    }
);

close $FILEHANDLE;

## File is created and has content
ok( -s $bash_file_path, q{Create bash} );

## Run command
my $ok = run( command => [ q{bash}, $bash_file_path ] );

## Created reformated output file
ok( -s $outfile_path, q{Created: } . $outfile_path );

## Create anonymous filehandle
$FILEHANDLE = IO::Handle->new();

## Open filehandle
open $FILEHANDLE, q{<}, $outfile_path
  or croak( q{Cannot read '} . $outfile_path . q{' :} . $OS_ERROR . $NEWLINE );

## Test the outfile from the bash script is properly formatted
_parse_outfile( { FILEHANDLE => $FILEHANDLE, } );

close $FILEHANDLE;

## Remove test files
unlink $bash_file_path, $outfile_path or croak q{Could not remove test files};

done_testing();

######################
####SubRoutines#######
######################

sub _build_test_file_recipe {

## Function : Builds the test file for testing the housekeeping function
## Returns  : ""
## Arguments: $bash_file_path             => Test file to write recipe to
##          : $commands_ref               => Commands to stream to perl oneliner
##          : $FILEHANDLE                 => Sbatch filehandle to write to
##          : $log_file_path              => The log file {REF}
##          : $reformat_sacct_headers_ref => Reformated sacct headers

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bash_file_path;
    my $commands_ref;
    my $FILEHANDLE;
    my $log_file_path;
    my $reformat_sacct_headers_ref;

    my $tmpl = {
        commands_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$commands_ref,
            strict_type => 1,
        },
        reformat_sacct_headers_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$reformat_sacct_headers_ref,
            strict_type => 1,
        },
        log_file_path => {
            store       => \$log_file_path,
            strict_type => 1,
        },
        FILEHANDLE => {
            required => 1,
            store    => \$FILEHANDLE,
        },
        bash_file_path => {
            required => 1,
            store    => \$bash_file_path,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw(Could not parse arguments!);

    # Add bash shebang
    build_shebang(
        {
            FILEHANDLE => $FILEHANDLE,
        }
    );

    ## Set shell attributes
    gnu_set(
        {
            FILEHANDLE  => $FILEHANDLE,
            set_errexit => 1,
            set_nounset => 1,
        }
    );

    slurm_reformat_sacct_output(
        {
            commands_ref               => \@commands,
            FILEHANDLE                 => $FILEHANDLE,
            log_file_path              => $log_file_path,
            reformat_sacct_headers_ref => \@reformat_sacct_headers,
        }
    );
    return;
}

sub _parse_outfile {

## Function : Test the outfile from the bash script is properly formatted
## Returns  : ""
## Aeguments: $FILEHANDLE => Filehandle to read from

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;

    my $tmpl = {
        FILEHANDLE => {
            required => 1,
            store    => \$FILEHANDLE,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    Readonly my $EXPECTED_ENTRIES => 11;
    Readonly my $JOB_ID           => 815_575;

    # Read file line by line
    while (<$FILEHANDLE>) {

        my $line = $_;

        chomp $line;

        #Header line
        if ( $NR == 1 ) {

            my @headers = split $TAB, $line;

            like( $line, qr/^#/xms, q{Found header line} );

            is( scalar @headers,
                $EXPECTED_ENTRIES, q{Checking number of expected headers} );

            is( $headers[0], q{#JobID}, q{Checking first header} );

            is( $headers[-1], q{ExitCode}, q{Checking last header} );
        }
        else {

            my @job_info_entries = split $TAB, $line;

            is( scalar @job_info_entries,
                $EXPECTED_ENTRIES, q{Checking number of expected entries} );

            is( $job_info_entries[0], $JOB_ID, q{Checking first entry} );

            is( $job_info_entries[-1], q{0:0}, q{Checking last entry} );
        }
    }
    return;
}
