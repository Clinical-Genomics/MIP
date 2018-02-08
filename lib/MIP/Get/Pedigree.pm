package MIP::Get::Pedigree;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use Readonly;
use List::MoreUtils qw{ any };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ get_sample_info };
}

sub get_sample_info {

    ## Function     : Given pedigree hash, extract an specific field given a sample level pedigree info
    ## Returns      : @info_out
    ## Arguments    : $get_values_for_key             =>  Key in pedigree hash to return its value
    ##              : $sample_info_intersect_key      =>  Key in pedigree hash to search
    ##              : $sample_info_intersect_value    =>  Value in pedigree hash for $sample_info_intersect_key
    ##              : $pedigree_hash                  =>  Pedigree hash {REF}
    ## Note         : It will search for "$sample_info_intersect_key: $sample_info_intersect_value" in pedigree_hash
    ##                and return the value for the key specified in $get_values_for_key.

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $get_values_for_key;
    my $pedigree_href;
    my $sample_info_intersect_key;
    my $sample_info_intersect_value;

    my $tmpl = {
        get_values_for_key => {
            defined     => 1,
            required    => 1,
            store       => \$get_values_for_key,
            strict_type => 1,
        },
        pedigree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$pedigree_href,
            strict_type => 1,
        },
        sample_info_intersect_key => {
            defined     => 1,
            store       => \$sample_info_intersect_key,
            strict_type => 1,
        },
        sample_info_intersect_value => {
            defined     => 1,
            store       => \$sample_info_intersect_value,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @info_out;

    if ( not $sample_info_intersect_key or not $sample_info_intersect_value ) {

        DEFAULT:
        foreach my $sample_href ( @{ $pedigree_href->{samples} } ) {

            ## Continue if get_values_for_key does not exist
            next DEFAULT if ( not $sample_href->{$get_values_for_key} );

            push @info_out, $sample_href->{$get_values_for_key};
        }
    }
    else {

      SAMPLE:
        foreach my $sample_href ( @{ $pedigree_href->{samples} } ) {

            ## Check if sample_info_intersect_key entry exists in this particular sample
            next SAMPLE if ( not $sample_href->{$sample_info_intersect_key} );

            ## Continue if sample_info_intersect_value is not matching
            next SAMPLE
              if (
                not $sample_href->{$sample_info_intersect_key} eq
                $sample_info_intersect_value );

            ## Continue if get_values_for_key does not exist
            next SAMPLE if ( not $sample_href->{$get_values_for_key} );

            if ( $sample_href->{$sample_info_intersect_key} eq
                $sample_info_intersect_value )
            {

                push @info_out, $sample_href->{$get_values_for_key};
            }

        }
    }

    return @info_out;
}

1;
