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

use FindBin qw($Bin);    #Find directory of script
use File::Basename qw(dirname basename);
use File::Spec::Functions qw(catdir);
use Getopt::Long;
use Test::More;

## Third party module(s)
use List::Util qw(any);

## MIPs lib/
use lib catdir( dirname($Bin), 'lib' );
use MIP::Script::Utils qw(help);
use MIP::Test::Writefile qw(test_write_to_file);

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = '0.0.0';

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
    help(
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

    $perl_module{'MIP::Script::Utils'} = [qw(help)];

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

use MIP::Workloadmanager::Slurm qw(slurm_build_sbatch_header);

my $NEWLINE = q{\n};

diag(
"Test slurm_build_sbatch_header $MIP::Workloadmanager::Slurm::VERSION, Perl $^V, $EXECUTABLE_NAME"
);

## Base arguments
my $sbatch_shebang = '#SBATCH ';

my %base_argument = (
    FILEHANDLE => {
        input           => undef,
        expected_output => $sbatch_shebang,
    },
);

## Specific arguments
my %argument = (
    project_id => {
        input           => 'project_id_test',
        expected_output => '--account=project_id_test',
    },
    core_number => {
        input           => '8',
        expected_output => '--ntasks=8',
    },
    process_time => {
        input           => '1:00:00',
        expected_output => '--time=1:00:00',
    },
    slurm_quality_of_service => {
        input           => 'high',
        expected_output => '--qos=high',
    },
    job_name => {
        input           => 'test_job_name',
        expected_output => '--job-name=test_job_name',
    },
    stdoutfile_path => {
        input           => 'outfile.test',
        expected_output => '--output=outfile.test',
    },
    stderrfile_path => {
        input           => 'stderrfile.test',
        expected_output => '--error=stderrfile.test',
    },
    email => {
        input           => 'test.testsson@test.com',
        expected_output => '--mail-user=test.testsson@test.com',
    },
    email_types_ref => {
        inputs_ref      => [ 'BEGIN', 'FAIL', 'END' ],
        expected_output => '--mail-type=BEGIN,FAIL,END',
    },
);

my @commands = slurm_build_sbatch_header(
    {
        project_id               => $argument{project_id}{input},
        slurm_quality_of_service => $argument{slurm_quality_of_service}{input},
        job_name                 => $argument{job_name}{input},
        stderrfile_path          => $argument{stderrfile_path}{input},
        stdoutfile_path          => $argument{stdoutfile_path}{input},
        email                    => $argument{email}{input},
        email_types_ref          => $argument{email_types_ref}{inputs_ref},
        process_time             => $argument{process_time}{input},
        core_number              => $argument{core_number}{input},
    }
);

## Testing return of commands
foreach my $key ( keys %argument ) {

    # Add sbatch shebang to all expected outputs
    my $expected_output = $sbatch_shebang . $argument{$key}{expected_output};

    ok( ( any { $_ eq $expected_output } @commands ), 'Argument: ' . $key );
}

## Testing write to file

# Fake arguments
my @args = (
    project_id => 1,
    FILEHANDLE => undef,
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&slurm_build_sbatch_header;

my $function_base_command = '#SBATCH ';

test_write_to_file(
    {
        args_ref             => \@args,
        module_function_cref => $module_function_cref,
        base_command         => $function_base_command,
        separator            => $NEWLINE,
    }
);

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
