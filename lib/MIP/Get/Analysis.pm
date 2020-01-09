package MIP::Get::Analysis;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## Third party module(s)
use autodie;
use List::MoreUtils qw{ all any };
use Log::Log4perl;
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $EMPTY_STR $NEWLINE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.15;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      get_overall_analysis_type
      get_vcf_parser_analysis_suffix
    };
}

sub get_overall_analysis_type {

## Function : Detect if all samples has the same sequencing type and return consensus or mixed
## Returns  : q{consensus} | q{mixed} - analysis_type
## Arguments: $analysis_type_href => Analysis_type hash {REF}
##          : $log                => Log Object

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $analysis_type_href;
    my $log;

    my $tmpl = {
        analysis_type_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$analysis_type_href,
            strict_type => 1,
        },
        log => {
            required => 1,
            store    => \$log,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @analysis_types = (qw{ dragen_rd_dna vrn wes wgs wts });

  ANALYSIS:
    foreach my $analysis_type (@analysis_types) {

        ## If consensus is reached
        if ( all { $_ eq $analysis_type } values %{$analysis_type_href} ) {

            return $analysis_type;
        }
    }

    ## Check that the user supplied analysis type is supported
    foreach my $user_analysis_type ( values %{$analysis_type_href} ) {

        if ( not any { $_ eq $user_analysis_type } @analysis_types ) {

            $log->fatal(
                q{'} . $user_analysis_type . q{' is not a supported analysis_type} );
            $log->fatal( q{Supported analysis types are '}
                  . join( q{', '}, @analysis_types )
                  . q(') );
            $log->fatal(q{Aborting run});
            exit 1;
        }
    }

    # No consensus, then it must be mixed
    return q{mixed};
}

sub get_vcf_parser_analysis_suffix {

## Function : Get the vcf parser analysis suffix
## Returns  : @analysis_suffixes
## Arguments: $vcfparser_outfile_count => Number of user supplied vcf parser outfiles

    my ($arg_href) = @_;

## Flatten argument(s)
    my $vcfparser_outfile_count;

    my $tmpl = {
        vcfparser_outfile_count => {
            defined     => 1,
            required    => 1,
            store       => \$vcfparser_outfile_count,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    Readonly my $VCFPARSER_OUTFILE_COUNT => $vcfparser_outfile_count - 1;

    my @analysis_suffixes;

    ## Determined by vcfparser output
    # Set research (="") and selected file suffix
    for my $vcfparser_outfile_counter ( 0 .. $VCFPARSER_OUTFILE_COUNT ) {

        if ( $vcfparser_outfile_counter == 1 ) {

            ## Select file variants
            push @analysis_suffixes, q{selected};
            next;
        }
        push @analysis_suffixes, $EMPTY_STR;
    }
    return @analysis_suffixes;
}

1;
