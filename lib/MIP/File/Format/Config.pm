package MIP::File::Format::Config;

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

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ write_mip_config };
}

sub write_mip_config {

## Function : Write config file for analysis
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $log                   => Log object
##          : $remove_keys_ref       => Keys to remove before writing to file {REF}
##          : $sample_info_href      => Info on samples and family hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $log;
    my $remove_keys_ref;
    my $sample_info_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        remove_keys_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$remove_keys_ref,
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

    use File::Basename qw{ dirname };
    use File::Path qw{ make_path };
    use MIP::File::Format::Yaml qw{ write_yaml };

    return if ( not $active_parameter_href->{config_file_analysis} );

    ## Create directory unless it already exists
    make_path( dirname( $active_parameter_href->{config_file_analysis} ) );

    ## Remove previous analysis specific info not relevant for current run e.g. log file, sample_ids which are read from pedigree or cmd
    delete @{$active_parameter_href}{ @{$remove_keys_ref} };

    ## Writes a YAML hash to file
    write_yaml(
        {
            yaml_href      => $active_parameter_href,
            yaml_file_path => $active_parameter_href->{config_file_analysis},
        }
    );
    $log->info( q{Wrote: } . $active_parameter_href->{config_file_analysis} );

    ## Add to sample_info for use downstream
    $sample_info_href->{config_file_analysis} =
      $active_parameter_href->{config_file_analysis};
    return;
}

1;
