package MIP::File::Format::Store;

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
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $LOG_NAME $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ define_analysis_files_to_store parse_store_files set_analysis_files_to_store };
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

    my %analysis_store_file = define_analysis_files_to_store(
        { active_parameter_href => $active_parameter_href, } );

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
    my @store_files = grep { !$seen{ $_->{path} }++ } ( reverse @{$store_files_ref} );
    $store_files_ref = [ reverse @store_files ];

    return $store_files_ref;
}

1;
