package MIP::Recipes::Install::Mip_scripts;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use List::Util qw{ none };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Path::Tiny qw{ path };
use utf8;
use warnings qw{ FATAL utf8 };
use warnings;

## CPAN
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $DOT $LOG_NAME $NEWLINE $SPACE $UNDERSCORE };
use MIP::Environment::Child_process qw{ child_process };
use MIP::Program::Gnu::Coreutils qw{ gnu_cp };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ install_mip_scripts };
}

sub install_mip_scripts {

## Function : Install mip_scripts
## Returns  :
## Arguments: $active_parameter_href => Active parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Install qw{ check_mip_executable };

    Readonly my $MOVE_DIRS_UP => 5;

    my $conda_environment_path = $active_parameter_href->{conda_environment_path};
    my @select_programs        = @{ $active_parameter_href->{select_programs} };

    return if ( none { $_ eq q{mip_scripts} } @select_programs );

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Get mip directory relative to this file since $Bin might have been set already
    my $mip_dir_path = path(__FILE__)->parent($MOVE_DIRS_UP);

    ## Define MIP scripts and yaml files
    my @mip_scripts = qw{ cpanfile mip };

    my %mip_sub_script = (
        utility_scripts => [qw{ calculate_af.pl max_af.pl }],
        t         => [qw{ mip_install.test mip_analyse_rd_dna.test mip_core.t mip_analysis.test }],
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
              qc_eval_metric_-v1.3-.yaml
              qc_regexp_-v1.26-.yaml
              rank_model_-v1.33-.ini
              svrank_model_-v1.8-.ini
              }
        ],
    );

    my @mip_directories = qw{ lib t definitions };

    $log->info(q{Installing MIP's perl scripts});

    ## Check if mip installation exists and is executable
    # mip is proxy for all mip scripts
    check_mip_executable(
        {
            conda_environment_path => $conda_environment_path,
        }
    );

    ## Create directories
  DIRECTORY:
    foreach my $directory ( keys %mip_sub_script ) {

        path( catdir( $conda_environment_path, q{bin}, $directory ) )->mkpath();
    }

  DIRECTORY:
    foreach my $directory (@mip_directories) {

        my @cp_cmds = gnu_cp(
            {
                force        => 1,
                infile_path  => catdir( $mip_dir_path,           $directory ),
                outfile_path => catdir( $conda_environment_path, q{bin} ),
                recursive    => 1,
            }
        );

        my %process_return = child_process(
            {
                commands_ref => \@cp_cmds,
                process_type => q{ipc_cmd_run},
            }
        );

        if ( not $process_return{success} ) {

            $log->fatal(q{Failed to copy mip_scripts});
            $log->logdie( $process_return{error_message} );
        }
    }

    ## Copy mip scripts and sub scripts to conda env and make executable
  SCRIPT:
    foreach my $script (@mip_scripts) {

        my $src_path = catfile( $mip_dir_path, $script );
        my $dst_path = catfile( $conda_environment_path, q{bin}, $script );
        path($src_path)->copy($dst_path);
        path($dst_path)->chmod(q{a+x});
    }

  DIRECTORY:
    foreach my $directory ( keys %mip_sub_script ) {

      SCRIPT:
        foreach my $script ( @{ $mip_sub_script{$directory} } ) {

            my $src_path = catfile( $mip_dir_path, $directory, $script );
            my $dst_path = catfile( $conda_environment_path, q{bin}, $directory, $script );

            path($src_path)->copy($dst_path);
            path($dst_path)->chmod(q{a+x});
        }
    }
    return 1;
}

1;
