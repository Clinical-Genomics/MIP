#!/usr/bin/env perl

use 5.018;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname  };
use File::Spec::Functions qw{ catdir };
use FindBin qw{ $Bin };
use Getopt::Long;
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use List::Util qw(any);
use Modern::Perl qw{ 2014 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Script::Utils qw{ help };
use MIP::Test::Writefile qw(test_write_to_file);

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = 1.0.0;

## Constants
Readonly my $COMMA   => q{,};
Readonly my $CORE_NR => 8;
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };

### User Options
GetOptions(

    # Display help text
    q{h|help} => sub {
        done_testing();
        say {*STDOUT} $USAGE;
        exit;
    },

    # Display version number
    q{v|version} => sub {
        done_testing();
        say {*STDOUT} $NEWLINE
          . basename($PROGRAM_NAME)
          . $SPACE
          . $VERSION
          . $NEWLINE;
        exit;
    },
    q{vb|verbose} => $VERBOSE,
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
## Modules with import
    my %perl_module = ( q{MIP::Script::Utils} => [qw{ help }], );

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

## Modules
    my @modules = (q{MIP::Workloadmanager::Slurm});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
  }

use MIP::Workloadmanager::Slurm qw{ slurm_build_sbatch_header };

my $separator = q{\n};

diag(   q{Test slurm_build_sbatch_header from SLURM.pm v}
      . $MIP::Workloadmanager::Slurm::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
. $EXECUTABLE_NAME );

## Base arguments
my $sbatch_shebang = q{#SBATCH };

my %base_argument = (
    FILEHANDLE => {
        input           => undef,
        expected_output => $sbatch_shebang,
    },
);

## Specific arguments
my %argument = (
    core_number => {
        input           => $CORE_NR,
        expected_output => q{--ntasks=} . $CORE_NR,
    },
    email => {
        input           => q{test.testsson@test.com},
        expected_output => q{--mail-user=test.testsson@test.com},
    },
    email_types_ref => {
        inputs_ref      => [qw{ BEGIN FAIL END }],
        expected_output => q{--mail-type=BEGIN,FAIL,END},
    },
    job_name => {
        input           => q{test_job_name},
        expected_output => q{--job-name=test_job_name},
    },
    process_time => {
        input           => q{1:00:00},
        expected_output => q{--time=1:00:00},
    },
    project_id => {
        input           => q{project_id_test},
        expected_output => q{--account=project_id_test},
    },
    slurm_quality_of_service => {
        input           => q{high},
        expected_output => q{--qos=high},
    },
    stderrfile_path => {
        input           => q{stderrfile.test},
        expected_output => q{--error=stderrfile.test},
    },
    stdoutfile_path => {
        input           => q{outfile.test},
        expected_output => q{--output=outfile.test},
    },
);

my @commands = slurm_build_sbatch_header(
    {
        core_number              => $argument{core_number}{input},
        email                    => $argument{email}{input},
        email_types_ref          => $argument{email_types_ref}{inputs_ref},
        job_name                 => $argument{job_name}{input},
        process_time             => $argument{process_time}{input},
        project_id               => $argument{project_id}{input},
        slurm_quality_of_service => $argument{slurm_quality_of_service}{input},
        stderrfile_path          => $argument{stderrfile_path}{input},
        stdoutfile_path          => $argument{stdoutfile_path}{input},
    }
);

## Testing return of commands
foreach my $key ( keys %argument ) {

    # Add sbatch shebang to all expected outputs
    my $expected_output = $sbatch_shebang . $argument{$key}{expected_output};

    ok( ( any { $_ eq $expected_output } @commands ), q{Argument: } . $key );
}

## Testing write to file

# Fake arguments
my @args = (
    project_id => 1,
    FILEHANDLE => undef,
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&slurm_build_sbatch_header;

my $function_base_command = q{#SBATCH };

test_write_to_file(
    {
        args_ref             => \@args,
        base_command         => $function_base_command,
        module_function_cref => $module_function_cref,
        separator            => $separator,
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
