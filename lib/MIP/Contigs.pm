package MIP::Contigs;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use autodie qw{ :all };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };

## MIPs lib/
use MIP::Constants qw{ %PRIMARY_CONTIG };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.04;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ set_contigs };

}

sub set_contigs {

## Function : Set contig prefix and contig names depending on reference used.
## Returns  :
## Arguments: $file_info_href => File info hash {REF}
##          : $version        => Version of the human genome reference

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $version;

    my $tmpl = {
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        version => {
            defined     => 1,
            required    => 1,
            store       => \$version,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File_info qw{ set_alt_loci_contigs set_bam_contigs set_primary_contigs };

    ## Make a modifiable copy for downstream use of global constant
    my %primary_contig_clone = Readonly::Clone %PRIMARY_CONTIG;

    ## Get hash of genome build version primary assembly contigs
    my %primary_contig;
    @primary_contig{ @{ $primary_contig_clone{$version}{contigs} } } = ();

    ## Set alternative loci contig set
    set_alt_loci_contigs(
        {
            alt_contig_set_name => q{alt_loci},
            file_info_href      => $file_info_href,
            primary_contig_href => \%primary_contig,
        }
    );

    ## Set contigs sets for primary assembly
    my @primary_contig_sets = qw{ contigs contigs_size_ordered };

  PRIMARY_CONTIG_SET:
    foreach my $contig_set (@primary_contig_sets) {

        set_primary_contigs(
            {
                file_info_href      => $file_info_href,
                primary_contigs_ref => \@{ $primary_contig_clone{$version}{$contig_set} },
                primary_contig_set_name => $contig_set,
            }
        );
    }

## Set contigs sets for bam level processing downstream
    my %bam_contig_set = (
        bam_contigs              => q{contigs},
        bam_contigs_size_ordered => q{contigs_size_ordered},
    );

  BAM_CONTIG_SET:
    while ( my ( $bam_contig_set, $contig_set ) = each %bam_contig_set ) {

        set_bam_contigs(
            {
                file_info_href      => $file_info_href,
                primary_contigs_ref => \@{ $primary_contig_clone{$version}{$contig_set} },
                bam_contig_set_name => $bam_contig_set,
            }
        );
    }
    return;
}

1;
