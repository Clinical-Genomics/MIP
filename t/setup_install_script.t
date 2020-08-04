#!/usr/bin/env perl

use 5.026;
use Carp;
use Cwd;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2018 };
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COLON $COMMA $DOT $SPACE };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.01;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Script::Setup_script} => [qw{ setup_install_script }],
        q{MIP::Test::Fixtures}       => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Script::Setup_script qw{ setup_install_script };

diag(   q{Test setup_install_script from Setup_script.pm v}
      . $MIP::Script::Setup_script::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Create test log
my $log = test_log( {} );

my %active_parameter = (
    bash_set_errexit         => 1,
    bash_set_nounset         => 1,
    core_number              => 1,
    email_types              => [qw{ FAIL }],
    process_time             => q{1:00:00},
    project_id               => q{cust000},
    slurm_quality_of_service => q{low},
);

### Given a request to create bash file
# Create anonymous filehandle
my $filehandle = IO::Handle->new();

# For storing info to write
my $file_content;

# Store file content in memory by using referenced variable
open $filehandle, q{>}, \$file_content
  or croak q{Cannot write to} . $SPACE . $file_content . $COLON . $SPACE . $OS_ERROR;

## When the sub is launched
trap {
    setup_install_script(
        {
            active_parameter_href => \%active_parameter,
            filehandle            => $filehandle,
            file_name             => q{test.sh},
            invoke_login_shell    => 0,
            log                   => $log,
            remove_dir            => catfile( cwd(), $DOT, q{test} ),
            sbatch_mode           => 0,
            set_errexit           => $active_parameter{bash_set_errexit},
            set_nounset           => $active_parameter{bash_set_nounset},
        }
    );
};

# Close the filehandle
close $filehandle;

## Then 'shebang' should  be part of the header
ok( $file_content =~ / \A [#]!\s+\/usr\/bin\/env\/bash /xms, q{Create bash file} );

### Given input to write bash file for sbatch submission
# Create anonymous filehandle
my $filehandle_sbatch = IO::Handle->new();

# For storing info to write
my $file_content_sbatch;

# Store file content in memory by using referenced variable
open $filehandle_sbatch, q{>}, \$file_content_sbatch
  or croak q{Cannot write to}
  . $SPACE
  . $file_content_sbatch
  . $COLON
  . $SPACE
  . $OS_ERROR;

## When the sub is launched
trap {
    setup_install_script(
        {
            active_parameter_href => \%active_parameter,
            filehandle            => $filehandle_sbatch,
            file_name             => q{test.sh},
            invoke_login_shell    => 1,
            log                   => $log,
            remove_dir            => catfile( cwd(), $DOT, q{test} ),
            sbatch_mode           => 1,
            set_errexit           => $active_parameter{bash_set_errexit},
            set_nounset           => $active_parameter{bash_set_nounset},
        }
    );
};

# Close the filehandle
close $filehandle_sbatch;

## Then '#SBATCH' should  be part of the header
ok( $file_content_sbatch =~ / ^[#]SBATCH /xms, q{Create sbatch headers} );

done_testing();
