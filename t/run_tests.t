#!/usr/bin/env perl

use Modern::Perl qw{ 2014 };
use warnings qw{ FATAL utf8 };
use autodie;
use 5.018;
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };

use FindBin qw{ $Bin };
use File::Basename qw{ dirname basename };
use File::Spec::Functions qw{ catdir catfile devnull };
use Getopt::Long;
use IPC::Cmd qw{ can_run run };

## CPANM
use Readonly;

##MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Check::Modules qw{ check_perl_modules };
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

BEGIN {

    require MIP::Check::Modules;

    ## Special case to initiate testing
    my @modules = (q{Test::More});

    # Evaluate that all modules required are installed
    check_perl_modules(
        {
            modules_ref  => \@modules,
            program_name => $PROGRAM_NAME,
        }
    );

    ##Initate tests
    say {*STDOUT} q{Initiate tests:};
    say {*STDOUT} qq{\n} . q{Testing perl modules and selected functions},
      qq{\n};

    ## More proper testing
    use Test::More;

    ##Modules with import
    my %perl_module;

    $perl_module{autodie}                  = [qw{ open close :all }];
    $perl_module{charnames}                = [qw{ :full :short }];
    $perl_module{Cwd}                      = [qw{ abs_path }];
    $perl_module{q{File::Basename}}        = [qw{ dirname basename }];
    $perl_module{q{File::Path}}            = [qw{ make_path remove_tree }];
    $perl_module{q{File::Spec::Functions}} = [qw{ catfile catdir devnull }];
    $perl_module{FindBin}                  = [qw{ $Bin }];
    $perl_module{q{List::Util}}            = [qw{ any all uniq }];
    $perl_module{q{IPC::Cmd}}              = [qw{ can_run run }];
    $perl_module{q{Modern::Perl}}          = [qw{ 2014 }];
    $perl_module{open}                     = [qw{ :encoding(UTF-8) :std }];
    $perl_module{q{Params::Check}}         = [qw{ check allow last_error }];
    $perl_module{warnings}                 = [qw{ FATAL utf8 }];

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {

        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load } . $module;
    }

    ##Modules
    @modules = qw{ Cwd Getopt::Long IO::Handle
      IPC::System::Simple Log::Log4perl
      Path::Iterator::Rule POSIX strict
      TAP::Harness Time::Piece utf8
      YAML warnings };

  MODULE:
    for my $module (@modules) {

        require_ok($module) or BAIL_OUT q{Cannot load } . $module;
    }
}

## Constants
Readonly my $COMMA   => q{,};
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };

my $config_file = catfile( dirname($Bin), qw(templates mip_config.yaml) );
my $VERBOSE = 1;
our $VERSION = '0.0.2';

###User Options
GetOptions(
    q{c|config_file:s} => \$config_file,

    # Display help text
    q{h|help} => sub { done_testing(); say {*STDOUT} $USAGE; exit; },

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

use TAP::Harness;
use Cwd;

## Test central perl modules and import functions
test_modules();

mip_scripts();

ok( can_run(q{prove}), q{Checking can run perl prove} );

diag(   q{Test run_tests.t version }
      . $VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

##Run tests files using run_test_files.txt manifest
my $cmds_ref =
  [ qw(prove - < ), catdir( dirname($Bin), qw(t config run_test_files.txt) ), ];
my ( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
  run( command => $cmds_ref, verbose => $VERBOSE );

## Test MIP execuation
$cmds_ref = [ qw(prove mip.t :: -c ), $config_file ];
( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
  run( command => $cmds_ref, verbose => $VERBOSE );

## Reached the end safely
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

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return <<"END_USAGE";
 $program_name [options]
    -c/--config_file YAML config file for analysis parameters (defaults to ../templates/mip_config.yaml")
    -vb/--verbose Verbose
    -h/--help Display this help message
    -v/--version Display version
END_USAGE
}

sub test_modules {

##test_modules

##Function : Test perl modules and functions
##Returns  :
##Arguments:
##         :

    use Cwd;
    ok( getcwd(), q{Cwd: Locate current working directory} );

    use FindBin qw{ $Bin };
    ok( defined $Bin, q{FindBin: Locate directory of script} );

    use File::Basename qw{ dirname };
    ok( dirname($Bin),
        q{File::Basename qw{ dirname }: Strip the last part of directory} );

    use File::Spec::Functions qw{ catfile catdir devnull };
    ok( catdir( dirname($Bin), q{t} ),
        q{File::Spec::Functions qw{ catdir }: Concatenate directories} );
    ok( catfile( $Bin, q{run_tests.t} ),
        q{File::Spec::Functions qw{ catfile }: Concatenate files} );
    ok(
        catfile( dirname( devnull() ), q{stdout} ),
        q{File::Spec::Functions qw{ devnull }: Use devnull}
    );

    use File::Basename qw{ basename };
    my $file_path = catfile( $Bin, q{run_tests.t} );
    ok( basename($file_path),
        q{File::Basename qw{ basename }: Strip directories} );

    use Cwd qw{ abs_path };
    ok(
        abs_path( catfile( $Bin, q{run_tests.pl} ) ),
        q{Cwd_abs_path: Add absolute path}
    );

    use File::Path qw{ make_path remove_tree };
    ok( make_path(q{TEST}), q{File::Path_make_path: Create path} );

    ##Clean-up
    ok( remove_tree(q{TEST}), q{File::Path_remove_tree: Remove path} );

    ##MIPs lib/
    use lib catdir( dirname($Bin), q{lib} );
    use MIP::File::Format::Yaml qw{ load_yaml };
    use YAML;
    my $yaml_file =
      catdir( dirname($Bin), qw{ templates 643594-miptest_pedigree.yaml } );
    ok( -f $yaml_file,
        q{YAML: File=} . $yaml_file . q{ in MIP/templates directory} );

    my $yaml = load_yaml( { yaml_file => $yaml_file, } );

    # Check that we got something
    ok( defined $yaml, q{YAML: Load File} );
    ok( Dump($yaml),   q{YAML: Dump file} );

    use Params::Check qw{ check allow last_error };
    use Log::Log4perl;
    ## Creates log
    my $log_file = catdir( dirname($Bin), qw{ templates mip_log.yaml } );
    ok( -f $log_file,
        q{Log::Log4perl: File=} . $log_file . q{ in MIP directory} );

    use MIP::Log::MIP_log4perl qw{ initiate_logger };
    ## Creates log object
    my $log = initiate_logger(
        {
            categories_ref => [qw{ TRACE ScreenApp }],
            file_path      => $log_file,
            log_name       => q{Run_tests},
        }
    );

    ok( $log->info(1),  q{Log::Log4perl: info} );
    ok( $log->warn(1),  q{Log::Log4perl: warn} );
    ok( $log->error(1), q{Log::Log4perl: error} );
    ok( $log->fatal(1), q{Log::Log4perl: fatal} );

    use Getopt::Long;
    push @ARGV, qw{ -verbose 2 };

    my $verbose = 1;
    ok(
        GetOptions( q{verbose:n} => \$verbose ),
        q{Getopt::Long: Get options call}
    );
    ok( $verbose == 2, q{Getopt::Long: Get options modified} );

    ## Check time
    use Time::Piece;
    my $date_time = localtime;
    ok( $date_time, q{localtime = } . $date_time );
    my $date_time_stamp = $date_time->datetime;
    ok( $date_time_stamp, q{datetime = } . $date_time );
    my $date = $date_time->ymd;
    ok( $date, q{ymd = } . $date );

    ## Locate name of script
    ok( $PROGRAM_NAME, q{Detect program name = } . $PROGRAM_NAME );

    ## Execution of programs
    use IPC::Cmd qw{ can_run run };
    ok( can_run(q{perl}), q{Can run IPC::Cmd} );
    ok(
        my $bool = IPC::Cmd->can_capture_buffer,
        q{IPC::Cmd can capture buffer}
    );

    return;
}

sub mip_scripts {

## Function : Test MIP kit completion
## Returns  :
## Arguments:

    my @mip_scripts = qw{ download_reference.pl
      mip_install.pl mip.pl qccollect.pl
      vcfparser.pl perl_install.pl };

  SCRIPT:
    foreach my $script (@mip_scripts) {

        is( -e catfile( dirname($Bin), $script ),
            1, q{Found MIP file: } . $script );
    }

    my %mip_sub_scripts = (
        utility_scripts =>
          [qw{ calculate_af.pl covplots_exome.R covplots_genome.R max_af.pl }],
        definitions =>
          [qw{ define_download_references.yaml define_parameters.yaml }],
        t         => [qw{ mip_install.t mip.t run_tests.t mip_analysis.t }],
        templates => [
            qw{ mip_config.yaml mip_travis_config.yaml
              643594-miptest_pedigree.yaml
              mip_log.yaml }
        ],
    );

    my @mip_directories =
      ( qw{ definitions templates lib }, catdir(qw{ t data }), );

  DIRECTORY:
    foreach my $directory (@mip_directories) {

        is( -e catfile( dirname($Bin), $directory ),
            1, q{Found MIP sub dir: } . $directory );
    }
  DIRECTORY:
    foreach my $directory ( keys %mip_sub_scripts ) {

      SCRIPT:
        foreach my $script ( @{ $mip_sub_scripts{$directory} } ) {

            is( -e catfile( dirname($Bin), $directory, $script ),
                1, q{Found MIP sub file: } . $script );
        }
    }
    return;
}
