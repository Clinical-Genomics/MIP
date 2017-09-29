package MIP::Set::File;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use autodie;
use Params::Check qw{ check allow last_error };

BEGIN {

    use base qw{ Exporter };
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ set_file_suffix set_merged_infile_prefix };
}

sub set_file_suffix {

## set_file_suffix

## Function : Set the current file suffix for this job id chain
## Returns  : "$file_suffix"
## Arguments: $parameter_href, $suffix_key, $job_id_chain, $file_suffix
##          : $parameter_href => Holds all parameters
##          : $suffix_key     => Suffix key
##          : $job_id_chain   => Job id chain for program
##          : $file_suffix    => File suffix

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $suffix_key;
    my $job_id_chain;
    my $file_suffix;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        suffix_key => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$suffix_key
        },
        job_id_chain => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$job_id_chain
        },
        file_suffix => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$file_suffix
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    $parameter_href->{$suffix_key}{$job_id_chain} = $file_suffix;

    return $file_suffix;
}

sub set_merged_infile_prefix {

## set_merged_infile_prefix

## Function : Set the merged infile prefix for sample id
## Returns  :
## Arguments: $file_info_href, $sample_id, $merged_infile_prefix
##          : $file_info_href       => File info hash {REF}
##          : $sample_id            => Sample id
##          : $merged_infile_prefix => Merged infile prefix

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $sample_id;
    my $merged_infile_prefix;

    my $tmpl = {
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href
        },
        sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id
        },
        merged_infile_prefix => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$merged_infile_prefix
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    $file_info_href->{$sample_id}{merged_infile} = $merged_infile_prefix;

    return;
}

1;
