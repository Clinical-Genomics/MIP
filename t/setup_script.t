#!/usr/bin/env perl

use 5.018;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp;
use FindBin qw{ $Bin };
use Getopt::Long;
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
use MIP::Log::MIP_log4perl qw{ initiate_logger };
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = '1.0.0';

## Constants
Readonly my $COMMA      => q{,};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

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
    my %perl_module = (
        q{MIP::Log::MIP_log4perl} => [qw{ initiate_logger }],
        q{MIP::Script::Utils}     => [qw{ help }],
    );

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }

## Modules
    my @modules = (q{MIP::Script::Setup_script});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
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

## Create temp logger
my $test_dir = File::Temp->newdir();
my $test_log_path = catfile( $test_dir, q{test.log} );

## Creates log object
my $log = initiate_logger(
    {
        file_path => $test_log_path,
        log_name  => q{TEST},
    }
);

# Create anonymous filehandle
my $FILEHANDLE = IO::Handle->new();

my $directory_id      = q{family_1};
my $call_type         = q{both};
my $test_program_name = q{bwa_mem};
my $temp_dir          = catfile($test_dir);

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
    slurm_quality_of_service        => q{low},
    source_environment_commands_ref => [qw{ source activate test }],
    temp_directory                  => $temp_dir,
);
my %job_id;

my ($file_path) = setup_script(
    {
        active_parameter_href => \%active_parameter,
        call_type             => $call_type,
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
  . $UNDERSCORE
  . $call_type
  . $DOT;
my $file_path_prefix = catfile( $program_script_directory_path,
    q{dry_run} . $UNDERSCORE . $file_name_prefix );
my $file_name_suffix   = $DOT . q{sh};
my $expected_file_path = $file_path_prefix . q{0} . $file_name_suffix;
is( $file_path, $expected_file_path, q{Genereted file path} );

close $FILEHANDLE;

done_testing();

######################
####SubRoutines#######
######################

sub build_usage {

## Function  : Build the USAGE instructions
## Returns   :
## Arguments : $program_name => Name of the script

    my ($arg_href) = @_;

    ## Default(s)
    my $program_name;

    my $tmpl = {
        program_name => {
            default     => basename($PROGRAM_NAME),
            store       => \$program_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return <<"END_USAGE";
 $program_name [options]
    -vb/--verbose Verbose
    -h/--help     Display this help message
    -v/--version  Display version
END_USAGE
}
