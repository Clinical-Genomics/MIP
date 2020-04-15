package MIP::Delete::List;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };

## CPANM
use autodie qw { :all };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.06;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ delete_male_contig delete_non_wes_contig };
}

sub delete_male_contig {

## Function : Delete contig chrY|Y from contigs array if no male or other found
## Returns  : @contigs
## Arguments: $contigs_ref      => Contigs array to update {REF}
##          : $contig_names_ref => Contig names to remove {REF}
##          : $found_male       => Male(s) was included in the analysis

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $contig_names_ref;
    my $contigs_ref;
    my $found_male;

    my $tmpl = {
        contigs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$contigs_ref,
            strict_type => 1,
        },
        contig_names_ref => {
            allow => [
                sub {
                    check_allowed_array_values(
                        {
                            allowed_values_ref => [qw{ Y chrY }],
                            values_ref         => $arg_href->{contig_names_ref},
                        }
                    );
                }
            ],
            default     => [qw{ Y }],
            store       => \$contig_names_ref,
            strict_type => 1,
        },
        found_male => {
            allow       => qr{\A \d+ \z}sxm,
            defined     => 1,
            required    => 1,
            store       => \$found_male,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Parameter qw{ check_allowed_array_values };
    use MIP::Contigs qw{ delete_contig_elements };

    my @contigs = @{$contigs_ref};

    ## Removes contigY|chrY from contigs if no males or 'other' found in analysis
    if ( not $found_male ) {

        @contigs = delete_contig_elements(
            {
                contigs_ref        => \@contigs,
                remove_contigs_ref => $contig_names_ref,
            }
        );
    }
    return @contigs;
}

sub delete_non_wes_contig {

## Function : Delete contig chrM|MT from contigs array if consensus analysis type is wes
## Returns  : @contigs
## Arguments: $analysis_type_href => Analysis_type hash {REF}
##          : $contigs_ref        => Contigs array to update {REF}
##          : $contig_names_ref   => Contig names to remove {REF}
##          : $log                => Log object

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $analysis_type_href;
    my $contigs_ref;
    my $contig_names_ref;
    my $log;

    my $tmpl = {
        analysis_type_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$analysis_type_href,
            strict_type => 1,
        },
        contigs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$contigs_ref,
            strict_type => 1,
        },
        contig_names_ref => {
            allow => [
                sub {
                    check_allowed_array_values(
                        {
                            allowed_values_ref => [qw{ M MT }],
                            values_ref         => $arg_href->{contig_names_ref},
                        }
                    );
                }
            ],
            default     => [qw{ M MT }],
            store       => \$contig_names_ref,
            strict_type => 1,
        },
        log => {
            required => 1,
            store    => \$log,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Analysis qw{ get_overall_analysis_type };
    use MIP::Check::Parameter qw{ check_allowed_array_values };
    use MIP::Contigs qw{ delete_contig_elements };

    ## Detect if all samples has the same sequencing type and return consensus if reached
    my $consensus_analysis_type = get_overall_analysis_type(
        {
            analysis_type_href => $analysis_type_href,
        }
    );

    my @contigs = @{$contigs_ref};

    ## Removes contigM|chrMT from contigs
    if ( $consensus_analysis_type eq q{wes} or $consensus_analysis_type eq q{mixed} ) {

        @contigs = delete_contig_elements(
            {
                contigs_ref        => \@contigs,
                remove_contigs_ref => $contig_names_ref,
            }
        );
    }
    return @contigs;
}

1;
