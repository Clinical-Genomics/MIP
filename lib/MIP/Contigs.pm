package MIP::Contigs;

use 5.026;
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
use List::Util qw{ any };

## MIPs lib/
use MIP::Constants qw{ $LOG_NAME %PRIMARY_CONTIG $SPACE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.05;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ check_select_file_contigs
      set_contigs
      sort_contigs_to_contig_set
    };

}

sub check_select_file_contigs {

## Function : Check that select file contigs is a subset of primary contigs
## Returns  :
## Arguments: $contigs_ref             => Primary contigs of the human genome reference
##          : $select_file_contigs_ref => Select file contigs

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $contigs_ref;
    my $select_file_contigs_ref;

    my $tmpl = {
        contigs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$contigs_ref,
            strict_type => 1,
        },
        select_file_contigs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$select_file_contigs_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Array::Utils qw{ array_minus };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Check that select file contigs are a subset of primary contigs
    my @unique_select_contigs =
      array_minus( @{$select_file_contigs_ref}, @{$contigs_ref} );

    if (@unique_select_contigs) {

        $log->fatal( q{Option 'vcfparser_select_file' contig(s): } . join $SPACE,
            @unique_select_contigs );
        $log->fatal(
            q{Is not a subset of the human genome reference contigs: } . join $SPACE,
            @{$contigs_ref} );
        exit 1;
    }
    return 1;
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

sub sort_contigs_to_contig_set {

## Function : Sorts array depending on reference array. NOTE: Only entries present in reference array will survive in sorted array.
## Returns  : @sorted_contigs
## Arguments: $consensus_analysis_type    => Consensus analysis_type {REF}
##          : $sort_contigs_ref           => Contigs to sort according to reference contig set
##          : $sort_reference_contigs_ref => Contigs to use as reference map when sorting

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $consensus_analysis_type;
    my $sort_contigs_ref;
    my $sort_reference_contigs_ref;

    my $tmpl = {
        consensus_analysis_type => {
            defined     => 1,
            required    => 1,
            store       => \$consensus_analysis_type,
            strict_type => 1,
        },
        sort_contigs_ref => {
            default  => [],
            required => 1,
            store    => \$sort_contigs_ref,
        },
        sort_reference_contigs_ref => {
            default  => [],
            defined  => 1,
            required => 1,
            store    => \$sort_reference_contigs_ref,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my @sorted_contigs;

    ## Sanity check
    if ( not @{$sort_contigs_ref} ) {

        $log->fatal(q{Nothing to sort in supplied sort contigs set });
        exit 1;
    }
    if ( not @{$sort_reference_contigs_ref} ) {

        $log->fatal(q{No contigs to use in supplied sort reference contig set });
        exit 1;
    }

    ## Sort the contigs depending on reference array
  REF_ELEMENT:
    foreach my $element ( @{$sort_reference_contigs_ref} ) {

        ## If present in hash of array to sort push to @sorted_contigs
        if ( any { $_ eq $element } @{$sort_contigs_ref} ) {

            push @sorted_contigs, $element;
        }
    }

    ## Test if all contigs collected from select file was sorted by reference contig array
    if ( @sorted_contigs
        and scalar @{$sort_contigs_ref} != scalar @sorted_contigs )
    {

      SORT_ELEMENT:
        foreach my $element ( @{$sort_contigs_ref} ) {

            ## If element is not part of array
            if ( not any { $_ eq $element } @sorted_contigs ) {

                ## Special case when analysing wes since Mitochondrial contigs have no baits in exome capture kits
                next SORT_ELEMENT
                  if (  $consensus_analysis_type eq q{wes}
                    and $element =~ / MT$ | M$ /sxm );

                $log->fatal( q{Could not detect 'contig'= }
                      . $element
                      . q{ from column 1 in '-vcfparser_select_file' in reference contigs collected from '-human_genome_reference'}
                );
                exit 1;
            }
        }
    }
    return @sorted_contigs;
}

1;
