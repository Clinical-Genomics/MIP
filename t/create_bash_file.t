#!/usr/bin/env perl

use 5.026;
use Carp;
use Cwd;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catfile catdir };
use File::Temp qw{ tempfile };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2014 };
use Readonly;
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $COLON => q{;};
Readonly my $COMMA => q{,};
Readonly my $DOT   => q{.};
Readonly my $SPACE => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Language::Shell} => [qw{ create_bash_file }],
        q{MIP::Test::Fixtures}  => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Language::Shell qw{ create_bash_file };

diag(   q{Test create_bash_file from Shell.pm v}
      . $MIP::Language::Shell::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Create test log
my $log = test_log( {} );

my %parameter = (
    slurm_quality_of_service => q{low},
    project_id               => q{cust000},
    core_number              => 1,
    email_types              => [qw{ FAIL }],
    process_time             => q{1:00:00},
    slurm_quality_of_service => q{ low },
    bash_set_errexit         => 1,
    bash_set_nounset         => 1,
);

### Given a request to create bash file
# Create anonymous filehandle
my $FILEHANDLE = IO::Handle->new();

# For storing info to write
my $file_content;

# Store file content in memory by using referenced variable
open $FILEHANDLE, q{>}, \$file_content
  or croak q{Cannot write to} . $SPACE . $file_content . $COLON . $SPACE . $OS_ERROR;

## When the sub is launched
trap {
    create_bash_file(
        {
            FILEHANDLE         => $FILEHANDLE,
            file_name          => q{test.sh},
            remove_dir         => catfile( cwd(), $DOT, q{test} ),
            invoke_login_shell => 0,
            log                => $log,
            parameter_href     => \%parameter,
            sbatch_mode        => 0,
            set_errexit        => $parameter{bash_set_errexit},
            set_nounset        => $parameter{bash_set_nounset},
        }
    );
};

# Close the filehandle
close $FILEHANDLE;

## Then '--account' should  be part of the header
ok( $file_content =~ / (\/usr\/bin\/env\/bash)$ /xms, q{Create bash file} );

### Given input to write bash file for sbatch submission
# Create anonymous filehandle
my $FILEHANDLE_SBATCH = IO::Handle->new();

# For storing info to write
my $file_content_sbatch;

# Store file content in memory by using referenced variable
open $FILEHANDLE_SBATCH, q{>}, \$file_content_sbatch
  or croak q{Cannot write to}
  . $SPACE
  . $file_content_sbatch
  . $COLON
  . $SPACE
  . $OS_ERROR;

## When the sub is launched
trap {
    create_bash_file(
        {
            FILEHANDLE         => $FILEHANDLE_SBATCH,
            file_name          => q{test.sh},
            remove_dir         => catfile( cwd(), $DOT, q{test} ),
            invoke_login_shell => 1,
            log                => $log,
            parameter_href     => \%parameter,
            sbatch_mode        => 1,
            set_errexit        => $parameter{bash_set_errexit},
            set_nounset        => $parameter{bash_set_nounset},
        }
    );
};

# Close the filehandle
close $FILEHANDLE_SBATCH;

## Then '--account' should  be part of the header
ok( $file_content_sbatch =~ / ^(\#SBATCH) /xms, q{Create sbatch headers} );

### Given faulty input when requesting sbatch headers
$parameter{project_id} = undef;

# Create anonymous filehandle
my $FILEHANDLE_FAULTY = IO::Handle->new();

# For storing info to write
my $file_content_faulty;

# Store file content in memory by using referenced variable
open $FILEHANDLE_FAULTY, q{>}, \$file_content_faulty
  or croak q{Cannot write to}
  . $SPACE
  . $file_content_faulty
  . $COLON
  . $SPACE
  . $OS_ERROR;

## When the sub is launched
trap {
    create_bash_file(
        {
            FILEHANDLE         => $FILEHANDLE_FAULTY,
            file_name          => q{test.sh},
            remove_dir         => catfile( cwd(), $DOT, q{test} ),
            invoke_login_shell => 1,
            log                => $log,
            parameter_href     => \%parameter,
            sbatch_mode        => 1,
            set_errexit        => $parameter{bash_set_errexit},
            set_nounset        => $parameter{bash_set_nounset},
        }
    );
};

# Close the filehandle
close $FILEHANDLE_FAULTY;

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit on missing project_id} );
like( $trap->stderr, qr/FATAL/xms, q{Throw fatal log message} );

done_testing();
