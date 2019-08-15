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
use IPC::Cmd qw(can_run run);
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

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Gnu::Bash}       => [qw{ gnu_set }],
        q{MIP::Gnu::Coreutils}  => [qw{ gnu_mkdir }],
        q{MIP::Language::Shell} => [qw{ build_shebang create_housekeeping_function }],
        q{MIP::Test::Fixtures}  => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Gnu::Bash qw{ gnu_set };
use MIP::Gnu::Coreutils qw{ gnu_mkdir };
use MIP::Language::Shell qw{ build_shebang create_housekeeping_function };

diag(   q{Test create_housekeeping_function from Shell.pm v}
      . $MIP::Language::Shell::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

# Create anonymous filehandle
my $FILEHANDLE = IO::Handle->new();

# Create housekeeping function test sbatch file
my $bash_file_path = catfile( cwd(), q{test_create_housekeeping_function.sh} );
my $log_file_path  = catdir( cwd(), q{test_create_housekeeping_function} );

# Temporary directory
my $temp_dir = catdir( cwd(), q{.test_create_housekeeping_function} );

# Open filehandle for bash file
open $FILEHANDLE, q{>}, $bash_file_path
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
        RECIPE_FILEHANDLE     => $FILEHANDLE,
        recipe_log_file_path  => $log_file_path,
        recipe_temp_dir       => $temp_dir,
    }
);
close $FILEHANDLE;
close $LOG_FH;

## Testing write to file
ok( -e $bash_file_path, q{Create bash} );

ok( can_run(q{bash}), q{Checking can run bash binary} );

my $cmds_ref = [ q{bash}, $bash_file_path ];
my ( $success, $error_message, $full_buf_ref, $stdout_buf_ref, $stderr_buf_ref ) =
  run( command => $cmds_ref, verbose => $VERBOSE );

## Testing housekeeping function
ok( !-d $temp_dir, q{Performed housekeeping} );

remove_tree( $bash_file_path, $log_file_path . q{.status} );

done_testing();

######################
####SubRoutines#######
######################

sub _build_test_file_recipe {

##Function : Builds the test file for testing the housekeeping function
##Returns  : ""
##Arguments: $recipe_bash_file_path => Test file to write recipe to
##         : $RECIPE_FILEHANDLE     => FILEHANDLE to write to
##         : $recipe_log_file_path  => Log file path
##         : $recipe_temp_dir       => Temporary directory to use for test

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $recipe_bash_file_path;
    my $RECIPE_FILEHANDLE;
    my $recipe_log_file_path;
    my $recipe_temp_dir;

    my $tmpl = {
        recipe_bash_file_path => { required => 1, store => \$recipe_bash_file_path, },
        RECIPE_FILEHANDLE     => { required => 1, store => \$RECIPE_FILEHANDLE, },
        recipe_log_file_path  => { required => 1, store => \$recipe_log_file_path, },
        recipe_temp_dir       => { required => 1, store => \$recipe_temp_dir, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Add bash shebang
    build_shebang(
        {
            FILEHANDLE => $RECIPE_FILEHANDLE,
        }
    );

    ## Set shell attributes
    gnu_set(
        {
            FILEHANDLE  => $RECIPE_FILEHANDLE,
            set_errexit => 1,
            set_nounset => 1,
        }
    );

    # Create dir to test removal later
    gnu_mkdir(
        {
            FILEHANDLE       => $RECIPE_FILEHANDLE,
            indirectory_path => $recipe_temp_dir,
            parents          => 1,
        }
    );
    say {$RECIPE_FILEHANDLE} $NEWLINE;

    # Create housekeeping fucntion to remove temp_dir
    create_housekeeping_function(
        {
            FILEHANDLE         => $RECIPE_FILEHANDLE,
            job_ids_ref        => [qw{job_id_test}],
            log_file_path      => $recipe_log_file_path,
            remove_dir         => $recipe_temp_dir,
            trap_function_name => q{finish},
        }
    );
    return;
}
