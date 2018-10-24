#!/usr/bin/env perl

use 5.018;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp;
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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = '1.0.1';

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $COMMA      => q{,};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Script::Setup_script} => [qw{ setup_script }],
        q{MIP::Test::Fixtures}       => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Script::Setup_script qw{ setup_script };

diag(   q{Test setup_script from Setup_script.pm v}
      . $MIP::Script::Setup_script::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log();

## Create temp logger
my $test_dir = File::Temp->newdir();

# Create anonymous filehandle
my $FILEHANDLE = IO::Handle->new();

my $directory_id      = q{family_1};
my $temp_dir          = catfile($test_dir);
my $test_program_name = q{bwa_mem};

my %active_parameter = (
    bash_set_errexit                => 1,
    bash_set_nounset                => 1,
    bash_set_pipefail               => 1,
    bwa_mem                         => 2,
    email_types                     => [qw{ BEGIN FAIL }],
    outaligner_dir                  => q{bwa},
    outdata_dir                     => catfile( $test_dir, q{test_data_dir} ),
    outscript_dir                   => catfile( $test_dir, q{test_script_dir} ),
    project_id                      => q{wamdu},
			sacct_format_fields => [qw{ jobid }],
    slurm_quality_of_service        => q{low},
    source_environment_commands_ref => [qw{ source activate test }],
			submission_profile => q{slurm},
    temp_directory                  => $temp_dir,
);
my %job_id;

my ($file_path) = setup_script(
    {
        active_parameter_href => \%active_parameter,
        directory_id          => $directory_id,
        email_types_ref       => [qw{ FAIL }],
        FILEHANDLE            => $FILEHANDLE,
        job_id_href           => \%job_id,
        log                   => $log,
        program_directory     => $active_parameter{outaligner_dir},
        program_name          => $test_program_name,
        sleep                 => 1,
        source_environment_commands_ref =>
          \@{ $active_parameter{source_environment_commands_ref} },
    }
);

my $program_script_directory_path = catdir( $active_parameter{outscript_dir},
    $directory_id, $active_parameter{outaligner_dir} );
my $file_name_prefix =
    $test_program_name
  . $UNDERSCORE
  . $directory_id
  . $DOT;
my $file_path_prefix = catfile( $program_script_directory_path,
    q{dry_run} . $UNDERSCORE . $file_name_prefix );
my $file_name_suffix   = $DOT . q{sh};
my $expected_file_path = $file_path_prefix . q{0} . $file_name_suffix;
is( $file_path, $expected_file_path, q{Generated file path} );

close $FILEHANDLE;

done_testing();
