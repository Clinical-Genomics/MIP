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
use Modern::Perl qw{ 2018 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.02;

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
use MIP::Program::Gnu::Bash qw{ gnu_set };
use MIP::Workloadmanager::Slurm qw{ slurm_reformat_sacct_output };

diag(   q{Test slurm_reformat_sacct_output from Slurm.pm v}
      . $MIP::Workloadmanager::Slurm::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

###  PREPROCESSING

## Initialize variables
# Test input from sacct command
my $test_data_file_path =
  catfile( $Bin, qw{ data test_data slurm_reformat_sacct_output_input.txt } );

# Header for output status file
my @reformat_sacct_headers =
  qw{ JobID JobName Account Partition AllocCPUS TotalCPU Elapsed Start End State ExitCode };

# Write test bash script to this file
my $bash_file_path = catfile( $Bin, q{test_slurm_reformat_sacct_output.sh} );

# File with sacct output to reformat
my $log_file_path = catfile( $Bin, q{slurm_reformat_sacct_output.txt} );

# Output file with reformated sacct output
my $outfile_path = $log_file_path . q{.status};

## Create test file
# Create anonymous filehandle
my $filehandle = IO::Handle->new();

# Open filehandle
open $filehandle, q{>}, $bash_file_path
  or croak( q{Cannot write to '} . $bash_file_path . q{' :} . $OS_ERROR . $NEWLINE );

# Write reformat command to bash file
my @commands = ( q{less}, $test_data_file_path );
_build_test_file_recipe(
    {
        bash_file_path             => $bash_file_path,
        commands_ref               => \@commands,
        filehandle                 => $filehandle,
        log_file_path              => $log_file_path,
        reformat_sacct_headers_ref => \@reformat_sacct_headers,
    }
);
close $filehandle;

# File is created and has content
ok( -s $bash_file_path, q{Create bash} );

### TESTING STARTS

## When the reformat code parses a given a test file with sacct output
# Reformat sacct output
my $ok = run( command => [ q{bash}, $bash_file_path ] );

# Created reformated output file
ok( -s $outfile_path, q{Created: } . $outfile_path );

# Create anonymous filehandle
$filehandle = IO::Handle->new();

# Open filehandle
open $filehandle, q{<}, $outfile_path
  or croak( q{Cannot read '} . $outfile_path . q{' :} . $OS_ERROR . $NEWLINE );

## Then the tests should be passed (check the _parse_outfile sub for more information on the tests)
_parse_outfile( { filehandle => $filehandle, } );

### CLEANUP
close $filehandle;

# Remove test files
unlink $bash_file_path, $outfile_path or croak q{Could not remove test files};

done_testing();

######################
####SubRoutines#######
######################

sub _build_test_file_recipe {

## Function : Builds the test file for testing the housekeeping function
## Returns  :
## Arguments: $bash_file_path             => Test file to write recipe to
##          : $commands_ref               => Commands to stream to perl oneliner
##          : $filehandle                 => Sbatch filehandle to write to
##          : $log_file_path              => The log file {REF}
##          : $reformat_sacct_headers_ref => Reformated sacct headers

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bash_file_path;
    my $commands_ref;
    my $filehandle;
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
        filehandle => {
            required => 1,
            store    => \$filehandle,
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
            filehandle => $filehandle,
        }
    );
    print {$filehandle} $NEWLINE;

    ## Set shell attributes
    gnu_set(
        {
            filehandle  => $filehandle,
            set_errexit => 1,
            set_nounset => 1,
        }
    );

    slurm_reformat_sacct_output(
        {
            commands_ref               => \@commands,
            filehandle                 => $filehandle,
            log_file_path              => $log_file_path,
            reformat_sacct_headers_ref => \@reformat_sacct_headers,
        }
    );
    return;
}

sub _parse_outfile {

## Function : Test the outfile from the bash script is properly formatted
## Returns  :
## Aeguments: $filehandle => Filehandle to read from

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;

    my $tmpl = {
        filehandle => {
            required => 1,
            store    => \$filehandle,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    Readonly my $EXPECTED_ENTRIES => 11;
    Readonly my $JOB_ID           => 815_575;

    # Read file line by line
    while (<$filehandle>) {

        my $line = $_;

        chomp $line;

        # Check header line
        if ( $NR == 1 ) {

            my @headers = split $TAB, $line;

            # Starts with "#"
            like( $line, qr/^#/xms, q{Found header line} );

            # The number of header entries are correct
            is( scalar @headers,
                $EXPECTED_ENTRIES, q{Checking number of expected headers} );

            # The first header entry says "#JobID"
            is( $headers[0], q{#JobID}, q{Checking first header} );

            # The last header entry says "ExitCode"
            is( $headers[-1], q{ExitCode}, q{Checking last header} );
        }

        # Check the job lines
        else {

            my @job_info_entries = split $TAB, $line;

            # The number of entries are correct
            is( scalar @job_info_entries,
                $EXPECTED_ENTRIES, q{Checking number of expected entries} );

            # The first entry is the expected job id
            is( $job_info_entries[0], $JOB_ID, q{Checking first entry} );

            # The last entry is the correct exit code for the test file
            is( $job_info_entries[-1], q{0:0}, q{Checking last entry} );
        }
    }
    return;
}
