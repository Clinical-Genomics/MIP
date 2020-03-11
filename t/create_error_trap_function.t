#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Path qw{ remove_tree };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use IPC::Cmd qw{ can_run run };
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
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 0;
our $VERSION = 1.00;

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

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Program::Gnu::Bash}      => [qw{ gnu_set }],
        q{MIP::Program::Gnu::Coreutils} => [qw{ gnu_mkdir gnu_rm }],
        q{MIP::Language::Shell} =>
          [qw{ build_shebang enable_trap create_error_trap_function }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Program::Gnu::Bash qw{ gnu_set };
use MIP::Program::Gnu::Coreutils qw{ gnu_mkdir gnu_rm };
use MIP::Language::Shell qw{ build_shebang enable_trap create_error_trap_function };

diag(   q{Test create_error_trap_function from Shell.pm v}
      . $MIP::Language::Shell::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

# Create anonymous filehandle
my $filehandle = IO::Handle->new();

# Shell create error trap file
my $bash_file_path = catfile( cwd(), q{test_create_error_trap_function.sh} );
my $log_file_path  = catdir( cwd(), q{test_create_error_trap_function} );

# Temporary directory
my $temp_dir = catdir( cwd(), qw(test_dir .test_create_error_trap_function) );

# Open filehandle
open $filehandle, q{>}, $bash_file_path
  or croak( q{Cannot write to '} . $bash_file_path . q{' :} . $OS_ERROR . $NEWLINE );

## Open filehandle for log file
open my $LOG_FH, q{>}, $log_file_path . q{.status}
  or croak( q{Cannot write to '} . $bash_file_path . q{' :} . $OS_ERROR . $NEWLINE );

# Touch file
say {$LOG_FH} q{Logging};

## Write to bash file
_build_test_file_recipe(
    {
        recipe_bash_file_path => $bash_file_path,
        recipe_filehandle     => $filehandle,
        recipe_log_file_path  => $log_file_path,
        recipe_temp_dir       => $temp_dir,
    }
);

close $filehandle;
close $LOG_FH;

## Testing write to file
ok( -e $bash_file_path, q{Create bash} );

ok( can_run(q{bash}), q{Checking can run bash binary} );

my $cmds_ref = [ q{bash}, $bash_file_path ];
my ( $success, $error_message, $full_buf_ref, $stdout_buf_ref, $stderr_buf_ref ) =
  run( command => $cmds_ref, verbose => $VERBOSE );

## Testing error trap function
ok( $stderr_buf_ref->[-1] =~ /[:] \s+ Unknown \s+ Error \s+ - \s+ ExitCode[=]/sxm,
    q{Performed error trap} );

remove_tree( $log_file_path . q{.status} );

done_testing();

######################
####SubRoutines#######
######################

sub _build_test_file_recipe {

## Function : Builds the test file for testing the housekeeping function
## Returns  :
## Arguments: $recipe_bash_file_path => Test file to write recipe to
##          : $recipe_filehandle     => filehandle to write to
##          : $recipe_log_file_path  => Log file path
##          : $recipe_temp_dir       => Temporary directory to use for test

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $recipe_bash_file_path;
    my $recipe_filehandle;
    my $recipe_log_file_path;
    my $recipe_temp_dir;

    my $tmpl = {
        recipe_bash_file_path => { required => 1, store => \$recipe_bash_file_path },
        recipe_filehandle     => { required => 1, store => \$recipe_filehandle },
        recipe_log_file_path  => { required => 1, store => \$recipe_log_file_path, },
        recipe_temp_dir       => { required => 1, store => \$recipe_temp_dir },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Add bash shebang
    build_shebang(
        {
            filehandle => $recipe_filehandle,
        }
    );

    ## Set shell attributes
    gnu_set(
        {
            filehandle  => $recipe_filehandle,
            set_errexit => 1,
            set_nounset => 1,
        }
    );

    enable_trap(
        {
            filehandle         => $recipe_filehandle,
            trap_function_call => q{previous_command="$BASH_COMMAND"},
            trap_signals_ref   => [qw{ DEBUG }],
        }
    );

    # Create housekeeping fucntion to remove temp_dir
    create_error_trap_function(
        {
            filehandle         => $recipe_filehandle,
            job_ids_ref        => [qw{job_id_test}],
            log_file_path      => $recipe_log_file_path,
            trap_function_name => q{error},
        }
    );

    # Remove batch file to make clean exit
    gnu_rm(
        {
            filehandle  => $recipe_filehandle,
            infile_path => $recipe_bash_file_path,
        }
    );

    # Create dir to test removal later
    gnu_mkdir(
        {
            filehandle       => $recipe_filehandle,
            indirectory_path => $recipe_temp_dir,
            parents          => 0,
        }
    );
    say {$recipe_filehandle} $NEWLINE;

    return;
}
