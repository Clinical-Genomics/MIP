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
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ get_sample_origin get_sample_id get_sample_info };
}

## Constants
Readonly my $COMMA   => q{,};
Readonly my $DOT     => q{.};
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };

sub get_sample_origin {

    ## Function     : Given pedigree hash and sample id, get pedigree sample origin.
    ## Returns      : $sample_origin
    ## Arguments    : $analysis_type    =>  Analysis type to retrieve information
    ##              : $pedigree_hash    =>  Pedigree hash {REF}
    ##              : $sample_id        =>  Sample id from pedigree hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $analysis_type;
    my $pedigree_href;
    my $sample_id;

    my $tmpl = {
        analysis_type => {
            defined     => 1,
            store       => \$analysis_type,
            strict_type => 1,
        },
        pedigree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$pedigree_href,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    foreach my $sample ( @{ $pedigree_href->{samples} } ) {

        ## Check if sample_origin entry exists in this particular sample
        next if ( not $sample->{sample_origin} );

        ## Check if analysis_type exists in this particular sample
        next if ( not $sample->{analysis_type} eq $analysis_type );

        ## Continue if sample_id is not matching
        next if ( not $sample->{sample_id} eq $sample_id );

        if ( $sample->{sample_id} eq $sample_id ) {
            my $sample_origin = $sample->{sample_origin};

            return $sample_origin;
        }

    }

    return $SPACE;
}

sub get_sample_id {

    ## Function     : Given pedigree hash and sample origin, get pedigree sample id.
    ## Returns      : @sample_id
    ## Arguments    : $analysis_type    =>  Analysis type to retrieve information
    ##              : $pedigree_hash    =>  Pedigree hash {REF}
    ##              : $sample_origin        =>  Sample origin from pedigree hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $analysis_type;
    my $pedigree_href;
    my $sample_origin;

    my $tmpl = {
        analysis_type => {
            defined     => 1,
            store       => \$analysis_type,
            strict_type => 1,
        },
        pedigree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$pedigree_href,
            strict_type => 1,
        },
        sample_origin => {
            defined     => 1,
            required    => 1,
            store       => \$sample_origin,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @sample_id;

    foreach my $sample ( @{ $pedigree_href->{samples} } ) {

        ## Check if sample_origin entry exists in this particular sample
        next if ( not $sample->{sample_origin} );

        ## Check if analysis_type exists in this particular sample
        next if ( not $sample->{analysis_type} eq $analysis_type );

        ## Continue if sample_origin is not matching
        next if ( not $sample->{sample_origin} eq $sample_origin );

        if ( $sample->{sample_origin} eq $sample_origin ) {
            push @sample_id, $sample->{sample_id};

        }

    }

    return @sample_id;
}

sub get_sample_info {

    ## Function     : Given pedigree hash, extract an specific field given a sample level pedigree info
    ## Returns      : @sample_info
    ## Arguments    : $pedigree_hash        =>  Pedigree hash {REF}
    ##              : $sample_info_key      =>  Key in pedigree hash to search
    ##              : $sample_info_value    =>  Key in pedigree hash to match
    ##              : $sample_info_out      =>  Key in pedigree hash to return

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_info_key;
    my $sample_info_out;
    my $sample_info_value;
    my $pedigree_href;

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

        if ( $sample->{$sample_info_key} eq $sample_info_value ) {

            push @info_out, $sample->{$sample_info_out};
        }

    }

    return @info_out;
}

1;
