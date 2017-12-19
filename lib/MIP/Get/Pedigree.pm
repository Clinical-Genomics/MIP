package MIP::Get::Pedigree;

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
    our @EXPORT_OK = qw{ get_pedigree_family_info };
}

## Constants
Readonly my $SPACE => q{ };

sub get_pedigree_family_info {

## Function : Get the pedigree family keys and values
## Returns  :
## Arguments: $pedigree_href    => YAML pedigree info hash {REF}
##          : $sample_info_href => Info on samples and family hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $pedigree_href;
    my $sample_info_href;

    my $tmpl = {
        pedigree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$pedigree_href,
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

    ### Add key and values for family level info
  KEY:
    foreach my $key ( keys %{$pedigree_href} ) {

        ## Do not add sample level info
        next KEY if ( $key eq q{samples} );

        $sample_info_href->{$key} = $pedigree_href->{$key};
    }
    return;
}

1;
