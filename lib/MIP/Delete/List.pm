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
use Readonly;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ delete_contig_elements  delete_male_contig delete_non_wes_contig };
}

## Constants
Readonly my $SPACE => q{ };

sub delete_contig_elements {

## Function : Removes contig elements from array and return new contig array while leaving orginal elements_ref untouched.
## Returns  : @cleansed_contigs
## Arguments: $elements_ref       => Array to remove an element from {REF}
##          : $remove_contigs_ref => Remove these contigs {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $elements_ref;
    my $remove_contigs_ref;

    my $tmpl = {
        elements_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$elements_ref,
            strict_type => 1,
        },
        remove_contigs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$remove_contigs_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Copy array for sequencial removal
    my @cleansed_contigs  = @{$elements_ref};
    my @contigs_to_remove = @{$remove_contigs_ref};

    ### Make sure that contig is removed independent of genome source i.e prefix or not
    ## If contigs has prefix
    if ( defined $elements_ref->[0]
        && $elements_ref->[0] =~ / ^chr /xsm )
    {

        ## And remove contigs has not
        if ( defined $remove_contigs_ref->[0]
            && $remove_contigs_ref->[0] !~ / ^chr /xsm )
        {

            ## Add prefix to remove contigs to match contigs prefix
            @contigs_to_remove = map { q{chr} . $_ } @contigs_to_remove;
        }
    }

  CONTIG:
    foreach my $remove_contig (@contigs_to_remove) {

        ## Keep order of contigs i.e. do not use hash look-up
        @cleansed_contigs = grep { $_ ne $remove_contig } @cleansed_contigs;

    }
    return @cleansed_contigs;
}

sub delete_male_contig {

## Function : Delete contig chrY|Y from contigs array if no male or other found
## Returns  : @contigs
## Arguments: $contigs_ref      => Contigs array to update {REF}
##          : $contig_names_ref => Contig names to remove {REF}
##          : $found_male       => Male was included in the analysis

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
            allow       => [ 0, 1 ],
            defined     => 1,
            required    => 1,
            store       => \$found_male,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Parameter qw{ check_allowed_array_values };

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

    use MIP::Check::Parameter qw{ check_allowed_array_values };
    use MIP::Get::Analysis qw{ get_overall_analysis_type };

    ## Detect if all samples has the same sequencing type and return consensus if reached
    my $consensus_analysis_type = get_overall_analysis_type(
        {
            analysis_type_href => $analysis_type_href,
            log                => $log,
        }
    );

    my @contigs = @{$contigs_ref};

    ## Removes contigM|chrMT from contigs
    if ( $consensus_analysis_type eq q{wes} or $consensus_analysis_type eq q{mixed} ) {

        @contigs = delete_contig_elements(
            {
                elements_ref       => \@contigs,
                remove_contigs_ref => $contig_names_ref,
            }
        );
    }
    return @contigs;
}

1;
