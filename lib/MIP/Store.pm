package MIP::Store;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

## MIPs lib/
use MIP::Constants qw{ $LOG_NAME };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ define_analysis_files_to_store
      define_qc_metrics_to_store
      parse_store_files
      set_analysis_files_to_store
      store_files
      store_metrics};
}

sub define_analysis_files_to_store {

## Function : Define analysis files to store from MAIN
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}

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

    my %analysis_store_file = (
        config => {
            format => q{meta},
            path   => $active_parameter_href->{config_file},
            tag    => q{config},
        },
        config_analysis => {
            format => q{meta},
            path   => $active_parameter_href->{config_file_analysis},
            tag    => q{config_analysis},
        },
        log => {
            format => q{meta},
            path   => $active_parameter_href->{log_file},
            tag    => q{log},
        },
        pedigree => {
            format => q{meta},
            path   => $active_parameter_href->{pedigree_file},
            tag    => q{pedigree},
        },
        pedigree_fam => {
            format => q{meta},
            path   => $active_parameter_href->{pedigree_fam_file},
            tag    => q{pedigree_fam},
        },
        references_info => {
            format => q{meta},
            path   => $active_parameter_href->{reference_info_file},
            tag    => q{references_info},
        },
        sample_info => {
            format => q{meta},
            path   => $active_parameter_href->{sample_info_file},
            tag    => q{sample_info},
        },
    );
    return %analysis_store_file;
}

sub define_qc_metrics_to_store {

## Function : Define qc metrics to store
## Returns  :
## Arguments:

    my ($arg_href) = @_;

    my %store_metrics = (
        AT_DROPOUT => {
            analysis_mode => q{sample},
            recipe_name   => q{collecthsmetrics},
        },
        FOLD_80_BASE_PENALTY => {
            analysis_mode => q{sample},
            recipe_name   => q{collecthsmetrics},
        },
        fraction_duplicates => {
            analysis_mode => q{sample},
            recipe_name   => q{markduplicates},
        },
        GC_DROPOUT => {
            analysis_mode => q{sample},
            recipe_name   => q{collecthsmetrics},
        },
        gender => {
            analysis_mode => q{sample},
            recipe_name   => q{chanjo_sexcheck},
        },
        MEAN_INSERT_SIZE => {
            analysis_mode => q{sample},
            recipe_name   => q{collectmultiplemetricsinsertsize},
        },
        MEAN_TARGET_COVERAGE => {
            analysis_mode => q{sample},
            recipe_name   => q{collecthsmetrics}
        },
        MEDIAN_TARGET_COVERAGE => {
            analysis_mode => q{sample},
            recipe_name   => q{collecthsmetrics}
        },
        PCT_INTERGENIC_BASES => {
            analysis_mode => q{sample},
            recipe_name   => q{collectrnaseqmetrics},
        },
        PCT_INTRONIC_BASES => {
            analysis_mode => q{sample},
            recipe_name   => q{collectrnaseqmetrics},
        },
        PCT_MRNA_BASES => {
            analysis_mode => q{sample},
            recipe_name   => q{collectrnaseqmetrics},
        },
        PCT_OFF_BAIT => {
            analysis_mode => q{sample},
            recipe_name   => q{collecthsmetrics},
        },
        PCT_TARGET_BASES_10X => {
            analysis_mode => q{sample},
            recipe_name   => q{collecthsmetrics},
        },
        PCT_TARGET_BASES_20X => {
            analysis_mode => q{sample},
            recipe_name   => q{collecthsmetrics},
        },
        PCT_TARGET_BASES_30x => {
            analysis_mode => q{sample},
            recipe_name   => q{collecthsmetrics},
        },
        PCT_TARGET_BASES_50x => {
            analysis_mode => q{sample},
            recipe_name   => q{collecthsmetrics},
        },
        PCT_TARGET_BASES_100x => {
            analysis_mode => q{sample},
            recipe_name   => q{collecthsmetrics},
        },
        percentage_mapped_reads => {
            analysis_mode => q{sample},
            recipe_name   => q{bamstats},
        },
        percentage_uniquely_mapped_reads => {
            analysis_mode => q{sample},
            recipe_name   => q{star_log},
        },
        raw_total_sequences => {
            analysis_mode => q{sample},
            recipe_name   => q{bamstats},
        },
        reads_mapped => {
            analysis_mode => q{sample},
            recipe_name   => q{bamstats},
        },
    );
    return %store_metrics;
}

sub parse_store_files {

## Function : Parse store files and remove old duplicates based on same path
## Returns  : $store_files_ref
## Arguments: $store_files_ref => Store files {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $store_files_ref;

    my $tmpl = {
        store_files_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$store_files_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Remove duplicates, keep most recent additions
    my %seen;
    my @store_files = grep { not $seen{ $_->{path} }++ } ( reverse @{$store_files_ref} );
    $store_files_ref = [ reverse @store_files ];

    return $store_files_ref;
}

sub set_analysis_files_to_store {

## Function : Set analysis files to store
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $sample_info_href      => Sample info hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $sample_info_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
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

    use MIP::Sample_info qw{ set_file_path_to_store };

    my %analysis_store_file =
      define_analysis_files_to_store( { active_parameter_href => $active_parameter_href, } );

  FILE:
    foreach my $file ( keys %analysis_store_file ) {

        set_file_path_to_store(
            {
                format           => $analysis_store_file{$file}{format},
                path             => $analysis_store_file{$file}{path},
                id               => $active_parameter_href->{case_id},
                recipe_name      => q{mip_analyse},
                sample_info_href => $sample_info_href,
                tag              => $analysis_store_file{$file}{tag},
            }
        );
    }
    return;
}

sub store_files {

## Function : Store files from the analysis
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $sample_info_href      => Sample info hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $sample_info_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
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

    use MIP::Io::Write qw{ write_to_file };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    set_analysis_files_to_store(
        {
            active_parameter_href => $active_parameter_href,
            sample_info_href      => $sample_info_href,
        }
    );

    ## Parse and write store array to file
    my %store_files = (
        files => parse_store_files(
            {
                store_files_ref => $sample_info_href->{files},
            }
        )
    );

    ## Writes a YAML hash to file
    write_to_file(
        {
            data_href => \%store_files,
            format    => q{yaml},
            path      => $active_parameter_href->{store_file},
        }
    );
    $log->info( q{Wrote: } . $active_parameter_href->{store_file} );

    return;
}

sub store_metrics {

## Function : Store metrics from the analysis to file
## Returns  :
## Arguments: $qc_data_href          => Metric data hash {REF}
##          : $sample_info_href      => Sample info hash {REF}
##          : $store_metrics_outfile => Path to write store metrics to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $qc_data_href;
    my $sample_info_href;
    my $store_metrics_outfile;

    my $tmpl = {
        qc_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_data_href,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        store_metrics_outfile => {
            store       => \$store_metrics_outfile,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Io::Write qw{ write_to_file };
    use MIP::Qc_data qw{ get_qc_metric };
    use MIP::Sample_info qw{ get_pedigree_sample_ids };

    return if ( not defined $store_metrics_outfile );

    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my %metrics_deliverable;
    my %store_metrics = define_qc_metrics_to_store();

    ## Unpack
    my @sample_ids = get_pedigree_sample_ids( { sample_info_href => $sample_info_href, } );

  METRIC:
    while ( my ( $metric_name, $metric_href ) = each %store_metrics ) {

        my @ids =
          $metric_href->{analysis_mode} eq q{sample} ? @sample_ids : ( $sample_info_href->{case} );

      ID:
        foreach my $id (@ids) {

            push @{ $metrics_deliverable{metrics} },
              get_qc_metric(
                {
                    header       => $metric_href->{header},
                    id           => $id,
                    input        => $metric_href->{input},
                    metric_name  => $metric_name,
                    metric_value => $metric_href->{metric_value},
                    qc_data_href => $qc_data_href,
                    recipe_name  => $metric_href->{recipe_name},
                }
              );
        }
    }

    ## Writes a YAML hash to file
    write_to_file(
        {
            data_href => \%metrics_deliverable,
            format    => q{yaml},
            path      => $store_metrics_outfile,
        }
    );
    $log->info( q{Wrote: } . $store_metrics_outfile );

    return @{ $metrics_deliverable{metrics} };
}
1;
