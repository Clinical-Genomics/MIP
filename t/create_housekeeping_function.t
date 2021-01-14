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
use IPC::Cmd qw{ can_run };
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
use MIP::Constants qw{ $COMMA $NEWLINE $SPACE };


my $VERBOSE = 0;


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Environment::Child_process} => [qw{ child_process }],
        q{MIP::Program::Gnu::Bash}         => [qw{ gnu_set }],
        q{MIP::Program::Gnu::Coreutils}    => [qw{ gnu_mkdir }],
        q{MIP::Language::Shell} => [qw{ build_shebang create_housekeeping_function }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Environment::Child_process qw{ child_process };
use MIP::Program::Gnu::Bash qw{ gnu_set };
use MIP::Program::Gnu::Coreutils qw{ gnu_mkdir };
use MIP::Language::Shell qw{ build_shebang create_housekeeping_function };

diag(   q{Test create_housekeeping_function from Shell.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a create housekeeping function test sbatch file
my $bash_file_path = catfile( cwd(), q{test_create_housekeeping_function.sh} );
my $log_file_path  = catdir( cwd(), q{test_create_housekeeping_function} );

my $temp_dir = catdir( cwd(), q{.test_create_housekeeping_function} );

open my $filehandle, q{>}, $bash_file_path
  or croak( q{Cannot write to '} . $bash_file_path . q{' :} . $OS_ERROR . $NEWLINE );

open my $LOG_FH, q{>}, $log_file_path . q{.status}
  or croak( q{Cannot write to '} . $bash_file_path . q{' :} . $OS_ERROR . $NEWLINE );

# Touch file
say {$LOG_FH} q{Logging};

## When building the recipe
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

## Then the recipe file should exist
ok( -e $bash_file_path, q{Create bash} );

## Then the recipe file should be able to be executed
ok( can_run(q{bash}), q{Checking can run bash binary} );

my $cmds_ref       = [ q{bash}, $bash_file_path ];
my %process_return = child_process(
    {
        commands_ref => $cmds_ref,
        process_type => q{ipc_cmd_run},
    }
);

## Then process should return success
ok( $process_return{success}, q{Executed bash recipe successfully} );

## Then the housekeeping function should have removed the temp dir
ok( !-d $temp_dir, q{Performed housekeeping} );

## Clean-up
remove_tree( $bash_file_path, $log_file_path . q{.status} );

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
        recipe_bash_file_path => { required => 1, store => \$recipe_bash_file_path, },
        recipe_filehandle     => { required => 1, store => \$recipe_filehandle, },
        recipe_log_file_path  => { required => 1, store => \$recipe_log_file_path, },
        recipe_temp_dir       => { required => 1, store => \$recipe_temp_dir, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    build_shebang(
        {
            filehandle => $recipe_filehandle,
        }
    );
    print {$filehandle} $NEWLINE;

    gnu_set(
        {
            filehandle  => $recipe_filehandle,
            set_errexit => 1,
            set_nounset => 1,
        }
    );

    # Create dir to test removal later
    gnu_mkdir(
        {
            filehandle       => $recipe_filehandle,
            indirectory_path => $recipe_temp_dir,
            parents          => 1,
        }
    );
    say {$recipe_filehandle} $NEWLINE;

    # Create housekeeping function to remove temp_dir
    create_housekeeping_function(
        {
            filehandle         => $recipe_filehandle,
            job_ids_ref        => [qw{ job_id_test }],
            log_file_path      => $recipe_log_file_path,
            remove_dir         => $recipe_temp_dir,
            trap_function_name => q{finish},
        }
    );
    return;
}
