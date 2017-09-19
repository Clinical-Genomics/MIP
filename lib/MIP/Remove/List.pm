package MIP::Remove::List;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
# Allow unicode characters in this script
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };

use Readonly;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ remove_contig_elements };
}

## Constants
Readonly my $SPACE => q{ };

sub remove_contig_elements {

##remove_contig_elements

##Function : Removes contig elements from array and return new contig array while leaving orginal elements_ref untouched.
##Returns  : @cleansed_contigs
##Arguments: $elements_ref, $remove_contigs_ref
##         : $elements_ref       => Array to remove an element from {REF}
##         : $remove_contigs_ref => Remove these contigs {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $elements_ref;
    my $remove_contigs_ref;

    my $tmpl = {
        elements_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$elements_ref
        },
        remove_contigs_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$remove_contigs_ref
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Copy array for sequencial removal
    my @cleansed_contigs = @{$elements_ref};
    my @contigs_to_remove = @{$remove_contigs_ref};

    if ($elements_ref->[0] =~/ ^chr /xsm) {

      # Make sure that contig is removed independent of genome source
      @contigs_to_remove = map { q{chr} . $_ } @contigs_to_remove;
    }

  CONTIGS:
    foreach my $remove_contig (@contigs_to_remove) {

      @cleansed_contigs = grep { $_ ne $remove_contig } @cleansed_contigs;
    }
    return @cleansed_contigs;
}

1;
