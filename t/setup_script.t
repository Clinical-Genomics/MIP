#!/usr/bin/env perl

use 5.026;
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
use Modern::Perl qw{ 2018 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $DOT $SPACE $UNDERSCORE };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.03;

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

Readonly my $TWO      => 2;
Readonly my $ULIMIT_N => 4096;

## Create temp logger
my $log = test_log( {} );

my $test_dir = File::Temp->newdir();

# Create anonymous filehandle
my $filehandle = IO::Handle->new();

## Given input to create script
my $directory_id     = q{case_1};
my $temp_dir         = catfile($test_dir);
my $test_recipe_name = q{bwa_mem};

my %active_parameter = (
    bash_set_errexit                => 1,
    bash_set_nounset                => 1,
    bash_set_pipefail               => 1,
    bwa_mem                         => $TWO,
    email_types                     => [qw{ BEGIN FAIL }],
    outdata_dir                     => catfile( $test_dir, q{test_data_dir} ),
    outscript_dir                   => catfile( $test_dir, q{test_script_dir} ),
    project_id                      => q{wamdu},
    sacct_format_fields             => [qw{ jobid }],
    slurm_quality_of_service        => q{low},
    source_environment_commands_ref => [qw{ conda activate test }],
    submission_profile              => q{slurm},
    temp_directory                  => $temp_dir,
    ulimit_n                        => $ULIMIT_N,
);
my %job_id;

my ($recipe_file_path) = setup_script(
    {
        active_parameter_href => \%active_parameter,
        directory_id          => $directory_id,
        email_types_ref       => [qw{ FAIL }],
        filehandle            => $filehandle,
        job_id_href           => \%job_id,
        recipe_directory      => $test_recipe_name,
        recipe_name           => $test_recipe_name,
        sleep                 => 1,
        source_environment_commands_ref =>
          \@{ $active_parameter{source_environment_commands_ref} },
    }
);

my $recipe_script_directory_path =
  catdir( $active_parameter{outscript_dir}, $directory_id, $test_recipe_name );
my $file_name_prefix = $test_recipe_name . $UNDERSCORE . $directory_id . $DOT;
my $file_path_prefix =
  catfile( $recipe_script_directory_path, q{dry_run} . $UNDERSCORE . $file_name_prefix );
my $file_name_suffix   = $DOT . q{sh};
my $expected_file_path = $file_path_prefix . q{0} . $file_name_suffix;

## Then the recipe file path should be returned
is( $recipe_file_path, $expected_file_path, q{Generated recipe file path} );

close $filehandle;

done_testing();
