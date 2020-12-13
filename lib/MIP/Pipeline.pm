package MIP::Pipeline;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

## MIPs lib/
use MIP::Constants qw{ $LOG_NAME $SPACE};

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ run_analyse_pipeline run_download_pipeline run_install_pipeline };
}

sub run_analyse_pipeline {

## Function : Run analysis pipeline recipe
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $broadcasts_ref          => Holds the parameters info for broadcasting later {REF}
##          : $consensus_analysis_type => Consensus analysis type
##          : $file_info_href          => File info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $order_parameters_ref    => Order of parameters (for structured output) {REF}
##          : $order_recipes_ref       => Order of recipes
##          : $parameter_href          => Parameter hash {REF}
##          : $sample_info_href        => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $broadcasts_ref;
    my $consensus_analysis_type;
    my $file_info_href;
    my $job_id_href;
    my $order_parameters_ref;
    my $order_recipes_ref;
    my $parameter_href;
    my $sample_info_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        broadcasts_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$broadcasts_ref,
            strict_type => 1,
        },
        consensus_analysis_type => {
            defined     => 1,
            required    => 1,
            store       => \$consensus_analysis_type,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        order_parameters_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$order_parameters_ref,
            strict_type => 1,
        },
        order_recipes_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$order_recipes_ref,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Recipes
    use MIP::Recipes::Pipeline::Analyse_dragen_rd_dna
      qw{ pipeline_analyse_dragen_rd_dna };
    use MIP::Recipes::Pipeline::Analyse_rd_dna qw{ pipeline_analyse_rd_dna };
    use MIP::Recipes::Pipeline::Analyse_rd_dna_panel qw{ pipeline_analyse_rd_dna_panel };
    use MIP::Recipes::Pipeline::Analyse_rd_rna qw{ pipeline_analyse_rd_rna };
    use MIP::Recipes::Pipeline::Analyse_rd_dna_vcf_rerun
      qw{ pipeline_analyse_rd_dna_vcf_rerun };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Create dispatch table of pipelines
    my %pipeline = (
        dragen_rd_dna => \&pipeline_analyse_dragen_rd_dna,
        mixed         => \&pipeline_analyse_rd_dna,
        panel         => \&pipeline_analyse_rd_dna_panel,
        vrn           => \&pipeline_analyse_rd_dna_vcf_rerun,
        wes           => \&pipeline_analyse_rd_dna,
        wgs           => \&pipeline_analyse_rd_dna,
        wts           => \&pipeline_analyse_rd_rna,
    );

    $log->info( q{Pipeline analysis type: } . $consensus_analysis_type );
    $pipeline{$consensus_analysis_type}->(
        {
            active_parameter_href => $active_parameter_href,
            broadcasts_ref        => $broadcasts_ref,
            file_info_href        => $file_info_href,
            job_id_href           => $job_id_href,
            log                   => $log,
            order_parameters_ref  => $order_parameters_ref,
            order_recipes_ref     => $order_recipes_ref,
            parameter_href        => $parameter_href,
            sample_info_href      => $sample_info_href,
        }
    );
    return;
}

sub run_download_pipeline {

## Function : Run download pipeline recipe
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this download hash {REF}

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

    use MIP::Recipes::Pipeline::Download qw{ pipeline_download };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    pipeline_download(
        {
            active_parameter_href => $active_parameter_href,
            temp_directory        => $active_parameter_href->{temp_directory},
        }
    );
    return;
}

sub run_install_pipeline {

## Function : Run install pipeline recipe
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this install hash {REF}

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

    use MIP::Recipes::Pipeline::Install qw{ pipeline_install };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    $log->info(
        q{Installing pipelines: } . join $SPACE,
        @{ $active_parameter_href->{pipelines} }
    );
    pipeline_install(
        {
            active_parameter_href => $active_parameter_href,
            quiet                 => $active_parameter_href->{quiet},
            verbose               => $active_parameter_href->{verbose},
        }
    );
    return;
}

1;
