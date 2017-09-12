#!/usr/bin/env perl

#### Copyright 2017 Henrik Stranneheim

use Modern::Perl qw(2014);
use warnings qw(FATAL utf8);
use autodie;
use 5.018;    #Require at least perl 5.18
use utf8;
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );
use Carp;
use English qw(-no_match_vars);
use Params::Check qw(check allow last_error);

use Cwd;
use FindBin qw($Bin);    #Find directory of script
use File::Basename qw(dirname basename);
use File::Spec::Functions qw(catdir catfile);
use Getopt::Long;
use Test::More;
use IPC::Cmd qw(can_run run);

## MIPs lib/
use lib catdir( dirname($Bin), 'lib' );
use Script::Utils qw(help);

our $USAGE = build_usage( {} );

my $VERBOSE = 0;
our $VERSION = '1.0.0';

###User Options
GetOptions(
    'h|help' => sub {
        done_testing();
        print {*STDOUT} $USAGE, "\n";
        exit;
    },    #Display help text
    'v|version' => sub {
        done_testing();
        print {*STDOUT} "\n" . basename($PROGRAM_NAME) . q{  } . $VERSION,
          "\n\n";
        exit;
    },    #Display version number
    'vb|verbose' => $VERBOSE,
  )
  or (
    done_testing(),
    Script::Utils::help(
        {
            USAGE     => $USAGE,
            exit_code => 1,
        }
    )
  );

BEGIN {

### Check all internal dependency modules and imports
##Modules with import
    my %perl_module;

    $perl_module{'Script::Utils'} = [qw(help)];

    while ( my ( $module, $module_import ) = each %perl_module ) {

        use_ok( $module, @{$module_import} )
          or BAIL_OUT 'Cannot load ' . $module;
    }

##Modules
    my @modules = ('MIP::Workloadmanager::Slurm');

    for my $module (@modules) {

        require_ok($module) or BAIL_OUT 'Cannot load ' . $module;
    }
}

use MIP::Language::Shell qw(build_shebang);
use MIP::Gnu::Bash qw(gnu_set);
use MIP::Workloadmanager::Slurm qw(slurm_reformat_sacct_output);
use MIP::Test::Commands qw(test_function);

diag(
"Test slurm_reformat_sacct_output $MIP::Workloadmanager::Slurm::VERSION, Perl $^V, $EXECUTABLE_NAME"
);

## Set up test data parameters
# Test input from sacct command
my $test_data_file_path = catfile(
    cwd(), qw(data test_data
      slurm_reformat_sacct_output_input.txt)
);

# Header for output
my @reformat_sacct_headers = qw(JobID JobName Account Partition
  AllocCPUS TotalCPU Elapsed Start
  End State ExitCode);

# Slurm reformat sacct output test file
my $log_file_path = 'slurm_reformat_sacct_output.txt';

my @commands = ( 'less', $test_data_file_path );

# Create anonymous filehandle
my $FILEHANDLE = IO::Handle->new();

# Slurm reformat sacct output test sbatch file
my $bash_file_path = catfile( cwd(), 'test_slurm_reformat_sacct_output.sh' );

# Open filehandle
open $FILEHANDLE, '>', $bash_file_path
  or
  croak( q{Cannot write to '} . $bash_file_path . q{' :} . $OS_ERROR . "\n" );

## Write to bash file
_build_test_file_recipe(
    {
        commands_ref               => \@commands,
        reformat_sacct_headers_ref => \@reformat_sacct_headers,
        FILEHANDLE                 => $FILEHANDLE,
        bash_file_path             => $bash_file_path,
        log_file_path              => $log_file_path,
    }
);

close $FILEHANDLE;

## Testing write to file
ok( -e $bash_file_path, 'Create bash' );

ok( can_run('bash'), 'Checking can run bash binary' );

my $cmds_ref = [ 'bash', $bash_file_path ];
my ( $success, $error_message, $full_buf_ref, $stdout_buf_ref, $stderr_buf_ref )
  = run( command => $cmds_ref, verbose => $VERBOSE );

my $outfile = $log_file_path . '.status';

ok( -e $log_file_path . '.status', q{Created: } . $log_file_path . '.status' );

# Create anonymous filehandle
$FILEHANDLE = IO::Handle->new();

open $FILEHANDLE, '<', $outfile
  or croak( q{Cannot read '} . $outfile . q{' :} . $OS_ERROR . "\n" );

## Test the outfile from the bash script is properly formatted
_parse_outfile( { FILEHANDLE => $FILEHANDLE, } );

close $FILEHANDLE;

ok( can_run('rm'), 'Checking can run rm binary' );

$cmds_ref = [ 'rm', $bash_file_path, $outfile ];
( $success, $error_message, $full_buf_ref, $stdout_buf_ref, $stderr_buf_ref ) =
  run( command => $cmds_ref, verbose => $VERBOSE );

done_testing();

######################
####SubRoutines#######
######################

sub build_usage {

##build_usage

##Function : Build the USAGE instructions
##Returns  : ""
##Arguments: $program_name
##         : $program_name => Name of the script

    my ($arg_href) = @_;

    ## Default(s)
    my $program_name;

    my $tmpl = {
        program_name => {
            default     => basename($PROGRAM_NAME),
            strict_type => 1,
            store       => \$program_name,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw(Could not parse arguments!);

    return <<"END_USAGE";
 $program_name [options]
    -vb/--verbose Verbose
    -h/--help Display this help message
    -v/--version Display version
END_USAGE
}

sub _build_test_file_recipe {

##_build_test_file_recipe

##Function : Builds the test file for testing the housekeeping function
##Returns  : ""
##Arguments: $commands_ref, reformat_sacct_headers_ref, $log_file_path, $FILEHANDLE, $bash_file_path
##         : $commands_ref               => Commands to stream to perl oneliner
##         : $reformat_sacct_headers_ref => Reformated sacct headers
##         : $log_file_path              => The log file {REF}
##         : $FILEHANDLE                 => Sbatch filehandle to write to
##         : $bash_file_path             => Test file to write recipe to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $commands_ref;
    my $reformat_sacct_headers_ref;
    my $FILEHANDLE;
    my $log_file_path;
    my $bash_file_path;

    my $tmpl = {
        commands_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$commands_ref
        },
        reformat_sacct_headers_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$reformat_sacct_headers_ref
        },
        log_file_path  => { strict_type => 1, store => \$log_file_path },
        FILEHANDLE     => { required    => 1, store => \$FILEHANDLE },
        bash_file_path => { required    => 1, store => \$bash_file_path },
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
            reformat_sacct_headers_ref => \@reformat_sacct_headers,
            log_file_path              => $log_file_path,
            FILEHANDLE                 => $FILEHANDLE,
        }
    );
    return;
}

sub _parse_outfile {

##_parse_outfile

##Function : Test the outfile from the bash script is properly formatted
##Returns  : ""
##Arguments: $FILEHANDLE
##         : $FILEHANDLE => Filehandle to read from

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;

    my $tmpl = { FILEHANDLE => { required => 1, store => \$FILEHANDLE }, };

    check( $tmpl, $arg_href, 1 ) or croak qw(Could not parse arguments!);

    # Read file line by line
    while (<$FILEHANDLE>) {

        my $line = $_;

        chomp $line;

        #Header line
        if ( $. == 1 ) {

            my @headers = split "\t", $line;

            ok( $line =~ /^#/, q{Found header line} );

            ok( scalar @headers == 11, q{Checking number of expected headers} );

            ok( $headers[0] eq '#JobID', q{Checking first header} );

            ok( $headers[-1] eq 'ExitCode', q{Checking last header} );
        }
        else {

            my @job_info_entries = split "\t", $line;

            ok(
                scalar @job_info_entries == 11,
                q{Checking number of expected entries}
            );

            ok( $job_info_entries[0] eq '815575', q{Checking first entry} );

            ok( $job_info_entries[-1] eq q{0:0}, q{Checking last entry} );
        }
    }
    return;
}
