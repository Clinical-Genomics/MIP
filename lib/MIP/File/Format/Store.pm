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
use MIP::Constants qw{ $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ define_analysis_files_to_store set_analysis_files_to_store };
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
            file_type => q{meta},
            path      => $active_parameter_href->{config_file},
        },
        config_analysis => {
            file_type => q{meta},
            path      => $active_parameter_href->{config_file_analysis},
        },
        log => {
            file_type => q{meta},
            path      => $active_parameter_href->{log_file},
        },
        pedigree => {
            file_type => q{meta},
            path      => $active_parameter_href->{pedigree_file},
        },
        pedigree_fam => {
            file_type => q{meta},
            path      => $active_parameter_href->{pedigree_fam_file},
        },
        references_info => {
            file_type => q{meta},
            path      => $active_parameter_href->{reference_info_file},
        },
        sample_info => {
            file_type => q{meta},
            path      => $active_parameter_href->{sample_info_file},
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

  FILE_TAG:
    foreach my $file_tag ( keys %analysis_store_file ) {

        set_file_path_to_store(
            {
                file_tag         => $file_tag,
                file_type        => $analysis_store_file{$file_tag}{file_type},
                path             => $analysis_store_file{$file_tag}{path},
                sample_info_href => $sample_info_href,
            }
        );
    }
    return;
}

1;
