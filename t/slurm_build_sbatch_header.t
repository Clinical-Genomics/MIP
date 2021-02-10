#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir };
use FindBin qw{ $Bin };
use List::Util qw{ any };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Modern::Perl qw{ 2018 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };

use MIP::Test::Writefile qw{ test_write_to_file };

## Constants
Readonly my $CORE_NR => 8;
Readonly my $GPU_NR  => 1;

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Program::Slurm}  => [qw{ slurm_build_sbatch_header }],
        q{MIP::Test::Writefile} => [qw{ test_write_to_file }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Slurm qw{ slurm_build_sbatch_header };

my $separator = q{\n};

diag(   q{Test slurm_build_sbatch_header from Slurm.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a sbatch shebang and parameters
## Base arguments
my $sbatch_shebang = q{#SBATCH} . $SPACE;

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
    gpu_number => {
        input           => $GPU_NR,
        expected_output => q{--gres=gpu:} . $GPU_NR,
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
        gpu_number               => $argument{gpu_number}{input},
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

    ## Then commands should match expected output
    ok( ( any { $_ eq $expected_output } @commands ), q{Argument: } . $key );
}

## Testing write to file

# Given mock arguments
my @args = (
    project_id => 1,
    filehandle => undef,
);

## Coderef - enables generalized use of generate call
my $module_function_cref = \&slurm_build_sbatch_header;

my @function_base_commands = ( q{#SBATCH}, );

test_write_to_file(
    {
        args_ref             => \@args,
        base_commands_ref    => \@function_base_commands,
        module_function_cref => $module_function_cref,
        separator            => $separator,
    }
);

done_testing();
