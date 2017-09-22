package MIP::Update::Contigs;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use FindBin qw{ $Bin };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ update_contigs_for_run };
}

sub update_contigs_for_run {

## update_contigs_for_run

## Function : Update contigs depending on settings in run
## Returns  :
## Arguments: $file_info_href, $analysis_type_href, $found_male
##          : $file_info_href     => File info hash {REF}
##          : $analysis_type_href => Analysis_type hash {REF}
##          : $found_male         => Male was included in the analysis

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $analysis_type_href;
    my $found_male;

    my $tmpl = {
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href
        },
        analysis_type_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$analysis_type_href
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

    use MIP::Delete::List qw{ delete_non_wes_contig delete_male_contig };

    ## Delete contig chrM|MT from contigs array if consensus analysis type is wes
    my @wes_contig_arrays = (
        \@{ $file_info_href->{contigs_size_ordered} },
        \@{ $file_info_href->{contigs} },
        \@{ $file_info_href->{select_file_contigs} },
    );

  ARRAY_REF:
    foreach my $array_ref (@wes_contig_arrays) {

        @{$array_ref} = delete_non_wes_contig(
            {
                analysis_type_href => $analysis_type_href,
                contigs_ref        => $array_ref,
            }
        );
    }

    my @male_contig_arrays = (
        \@{ $file_info_href->{contigs_sv_size_ordered} },
        \@{ $file_info_href->{contigs_sv} },
        \@{ $file_info_href->{select_file_contigs} },
    );

  ARRAY_REF:
    foreach my $array_ref (@male_contig_arrays) {

        ## Removes contig_names from contigs array if no male or 'other' found
        @{$array_ref} = delete_male_contig(
            {
                contigs_ref => $array_ref,
                found_male  => $found_male,
            }
        );
    }
    return;
}

1;
