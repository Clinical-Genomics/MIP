package MIP::Delete::List;

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
    our @EXPORT_OK = qw{ delete_contig_elements delete_male_contig};
}

## Constants
Readonly my $SPACE => q{ };

sub delete_male_contig {

##delete_male_contig

##Function : Delete contig chrY|Y from contigs array if no male or other found
##Returns  : @contigs
##Arguments: $contigs_ref, $contig_names_ref, found_male
##         : $contigs_ref      => Contigs array to update {REF}
##         : $contig_names_ref => Contig names to remove {REF}
##         : $found_male       => Male was included in the analysis

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $contigs_ref;
    my $contig_names_ref;
    my $found_male;

    use MIP::Check::Parameter qw{ check_allowed_array_values };

    my $tmpl = {
        contigs_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$contigs_ref
        },
        contig_names_ref => {
            default => [qw{ Y }],
            allow   => [
                sub {
                    check_allowed_array_values(
                        {
                            allowed_values_ref => [qw{ Y chrY }],
                            values_ref         => $arg_href->{contig_names_ref},
                        }
                    );
                }
            ],
            strict_type => 1,
            store       => \$contig_names_ref
        },
        found_male => {
            required    => 1,
            defined     => 1,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$found_male,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @contigs = @{$contigs_ref};
    ## Removes contigY|chrY from contigs if no males or 'other' found in analysis
    if ( not $found_male ) {

        @contigs = delete_contig_elements(
            {
                elements_ref       => \@contigs,
                remove_contigs_ref => $contig_names_ref,
            }
        );
    }
    return @contigs;
}

sub delete_contig_elements {

##delete_contig_elements

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
    my @contigs_to_remove = @{$remove_contigs_ref};

    if ( $elements_ref->[0] =~ / ^chr /xsm ) {

        # Make sure that contig is removed independent of genome source
        @contigs_to_remove = map { q{chr} . $_ } @contigs_to_remove;
    }

    my %contig;

    # Remake allowed values into keys in contig hash
    map { $contig{$_} = undef } @{$elements_ref};

  CONTIGS:
    foreach my $remove_contig (@contigs_to_remove) {

        if ( exists $contig{$remove_contig} ) {

            delete $contig{$remove_contig};
        }
    }
    my @cleansed_contigs = keys %contig;

    return @cleansed_contigs;
}

1;
