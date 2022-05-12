#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use Getopt::Long;
use File::Basename qw{ dirname basename };
use File::Spec::Functions qw{ catdir catfile devnull };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use warnings qw{ FATAL utf8 };
use utf8;

## CPANM
use autodie;
use IPC::Cmd qw{ can_run run };
use Modern::Perl qw{ 2018 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $MIP_VERSION $NEWLINE $SPACE };
use MIP::Language::Perl qw{ check_modules_existance };
use MIP::Script::Utils qw{ help };

my $VERBOSE = 1;
our $USAGE = build_usage( {} );

BEGIN {

    require MIP::Language::Perl;

    ## Special case to initiate testing
    my @modules = (q{Test::More});

    # Evaluate that all modules required are installed
    check_modules_existance(
        {
            modules_ref  => \@modules,
            program_name => $PROGRAM_NAME,
        }
    );

    ##Initate tests
    say {*STDOUT} q{Initiate tests:};
    say {*STDOUT} qq{\n} . q{Testing perl modules and selected functions}, qq{\n};

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

    ## Modules
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

my $config_file = catfile( dirname($Bin), qw(templates mip_rd_dna_config.yaml) );

###User Options
GetOptions(
    q{c|config_file:s} => \$config_file,

    # Display help text
    q{h|help} => sub { done_testing(); say {*STDOUT} $USAGE; exit; },

    # Display version number
    q{v|version} => sub {
        done_testing();
        say {*STDOUT} $NEWLINE . basename($PROGRAM_NAME) . $SPACE . $MIP_VERSION . $NEWLINE;
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

diag(   q{Test mip_core.t version }
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Test central perl modules and import functions
test_modules();

mip_scripts();

ok( can_run(q{prove}), q{Checking can run perl prove} );

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
    -c/--config_file YAML config file for analysis parameters (defaults to ../templates/mip_rd_dna_config.yaml")
    -vb/--verbose Verbose
    -h/--help Display this help message
    -v/--version Display version
END_USAGE
}

sub test_modules {

##Function : Test perl modules and functions
##Returns  :
##Arguments:
##         :

    use Cwd;
    ok( getcwd(), q{Cwd: Locate current working directory} );

    use FindBin qw{ $Bin };
    ok( defined $Bin, q{FindBin: Locate directory of script} );

    use File::Basename qw{ dirname };
    ok( dirname($Bin), q{File::Basename qw{ dirname }: Strip the last part of directory} );

    use File::Spec::Functions qw{ catfile catdir devnull };
    ok( catdir( dirname($Bin), q{t} ),
        q{File::Spec::Functions qw{ catdir }: Concatenate directories} );
    ok( catfile( $Bin, q{mip_core.t} ), q{File::Spec::Functions qw{ catfile }: Concatenate files} );
    ok(
        catfile( dirname( devnull() ), q{stdout} ),
        q{File::Spec::Functions qw{ devnull }: Use devnull}
    );

    use File::Basename qw{ basename };
    my $file_path = catfile( $Bin, q{mip_core.t} );
    ok( basename($file_path), q{File::Basename qw{ basename }: Strip directories} );

    use Cwd qw{ abs_path };
    ok( abs_path( catfile( $Bin, q{mip_core.pl} ) ), q{Cwd_abs_path: Add absolute path} );

    use File::Path qw{ make_path remove_tree };
    ok( make_path(q{TEST}), q{File::Path_make_path: Create path} );

    ##Clean-up
    ok( remove_tree(q{TEST}), q{File::Path_remove_tree: Remove path} );

    ##MIPs lib/
    use lib catdir( dirname($Bin), q{lib} );
    use MIP::Io::Read qw{ read_from_file };
    use YAML;
    my $yaml_file = catdir( dirname($Bin), qw{ templates 643594-miptest_pedigree.yaml } );
    ok( -f $yaml_file, q{YAML: File=} . $yaml_file . q{ in MIP/templates directory} );

    my $yaml = read_from_file(
        {
            format => q{yaml},
            path   => $yaml_file,
        }
    );

    # Check that we got something
    ok( defined $yaml, q{YAML: Load File} );
    ok( Dump($yaml),   q{YAML: Dump file} );

    use Params::Check qw{ check allow last_error };
    use Log::Log4perl;
    ## Creates log
    my $log_file = catdir( dirname($Bin), qw{ templates mip_log.yaml } );
    ok( -f $log_file, q{Log::Log4perl: File=} . $log_file . q{ in MIP directory} );

    use MIP::Log::MIP_log4perl qw{ initiate_logger };
    ## Creates log object
    my $log = initiate_logger(
        {
            categories_ref => [qw{ TRACE ScreenApp }],
            file_path      => $log_file,
            log_name       => q{Mip_core},
        }
    );

    ok( $log->info(1),  q{Log::Log4perl: info} );
    ok( $log->warn(1),  q{Log::Log4perl: warn} );
    ok( $log->error(1), q{Log::Log4perl: error} );
    ok( $log->fatal(1), q{Log::Log4perl: fatal} );

    use Getopt::Long;
    push @ARGV, qw{ -verbose 2 };

    my $verbose = 1;
    ok( GetOptions( q{verbose:n} => \$verbose ), q{Getopt::Long: Get options call} );
    ok( $verbose == 2,                           q{Getopt::Long: Get options modified} );

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
    ok( can_run(q{perl}),                        q{Can run IPC::Cmd} );
    ok( my $bool = IPC::Cmd->can_capture_buffer, q{IPC::Cmd can capture buffer} );

    return;
}

sub mip_scripts {

## Function : Test MIP kit completion
## Returns  :
## Arguments:

    my @mip_scripts = qw{ mip cpanfile };

  SCRIPT:
    foreach my $script (@mip_scripts) {

        is( -e catfile( dirname($Bin), $script ), 1, q{Found MIP file: } . $script );
    }

    my %mip_sub_scripts = (
        utility_scripts => [qw{ calculate_af.pl max_af.pl }],
        definitions     => [
            qw{ analyse_parameters.yaml
              download_parameters.yaml
              install_parameters.yaml
              required_parameters.yaml
              mip_parameters.yaml
              not_required_parameters.yaml
              rd_dna_initiation_map.yaml
              rd_dna_parameters.yaml
              rd_rna_parameters.yaml
              rd_rna_initiation_map.yaml
              rd_dna_vcf_rerun_initiation_map.yaml
              rd_dna_vcf_rerun_parameters.yaml
              }
        ],
        t => [
            qw{ mip_install.test
              mip_analyse_rd_dna.test
              mip_analyse_rd_rna.test
              mip_analyse_rd_dna_vcf_rerun.test
              mip_core.t
              mip_analysis.test
              }
        ],
        templates => [
            qw{ 643594-miptest_pedigree.yaml
              gene_panels.bed
              grch38_mip_rd_dna_config.yaml
              mip_download_rd_dna_config_-1.0-.yaml
              mip_download_rd_rna_config_-1.0-.yaml
              mip_dragen_rd_dna_config.yaml
              mip_install_config.yaml
              mip_log.yaml
              mip_rd_dna_config.yaml
              mip_rd_dna_vcf_rerun_config.yaml
              mip_rd_rna_config.yaml
              program_test_cmds.yaml
              qc_regexp_-v1.26-.yaml
              rank_model_-v1.33-.ini
              svrank_model_-v1.8-.ini
              }
        ],
    );

    my @mip_directories = ( qw{ definitions templates lib }, catdir(qw{ t data }), );

  DIRECTORY:
    foreach my $directory (@mip_directories) {

        is( -e catfile( dirname($Bin), $directory ), 1, q{Found MIP sub dir: } . $directory );
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
