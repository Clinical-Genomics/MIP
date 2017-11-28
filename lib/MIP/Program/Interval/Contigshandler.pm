package MIP::Program::Interval::Contigshandler;

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
use Readonly;

## MIPs lib/
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ remove_array_element };
}

## Constants
Readonly my $SPACE => q{ };

sub remove_array_element {

## Function : Removes contigs from supplied contigs_ref.
## Returns  :
## Arguments: $contigs_ref            => The select file contigs {REF}
##          : $remove_contigs_ref     => Contig to be removed

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $contigs_ref;
    my $remove_contigs_ref;

    my $tmpl = {
        contigs_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$contigs_ref
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

    for ( my $index = 0 ; $index < scalar @{$contigs_ref} ; $index++ ) {

        foreach my $remove_contig ( @{$remove_contigs_ref} ) {

            if (   ( $contigs_ref->[$index] eq $remove_contig )
                || ( $contigs_ref->[$index] eq q{chr} . $remove_contig ) )
            {

                #Remove $element from array
                splice @{$contigs_ref}, $index, 1 ;
            }
        }
    }

    return;
}

1;
