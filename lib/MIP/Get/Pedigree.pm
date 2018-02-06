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
    ## Returns      : @sample_info
    ## Arguments    : $pedigree_hash        =>  Pedigree hash {REF}
    ##              : $sample_info_key      =>  Key in pedigree hash to search
    ##              : $sample_info_value    =>  Key in pedigree hash to match
    ##              : $sample_info_out      =>  Key in pedigree hash to return
    ## Note         : Returns an empty array if a particular key is not found.

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $pedigree_href;
    my $sample_info_key;
    my $sample_info_out;
    my $sample_info_value;

    my $tmpl = {
        pedigree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$pedigree_href,
            strict_type => 1,
        },
        sample_info_key => {
            defined     => 1,
            required    => 1,
            store       => \$sample_info_key,
            strict_type => 1,
        },
        sample_info_out => {
            defined     => 1,
            required    => 1,
            store       => \$sample_info_out,
            strict_type => 1,
        },
        sample_info_value => {
            defined     => 1,
            required    => 1,
            store       => \$sample_info_value,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @info_out;

    foreach my $sample ( @{ $pedigree_href->{samples} } ) {

        ## Check if sample_info_key entry exists in this particular sample
        next if ( not $sample->{$sample_info_key} );

        ## Continue if sample_info_value is not matching
        next if ( not $sample->{$sample_info_key} eq $sample_info_value );

        ## Continue if sample_info_out does not exist
        next if ( not $sample->{$sample_info_out});

        if ( $sample->{$sample_info_key} eq $sample_info_value ) {

            push @info_out, $sample->{$sample_info_out};
        }

    }

    return @info_out;
}

1;
