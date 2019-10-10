package MIP::Parse::Reference;

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
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $SPACE $TAB };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ parse_reference_for_vt };
}

sub parse_reference_for_vt {

## Function : Parse reference to make sure that they have been decomposed and normalised
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $log                     => Log object to write to
##          : $parameter_href          => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $log;
    my $parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        infile_lane_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_lane_prefix_href,
            strict_type => 1,
        },
        job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Reference qw{ check_references_for_vt };
    use MIP::Recipes::Analysis::Vt_core qw{ analysis_vt_core };

    return if ( not $active_parameter_href->{vt_decompose} );

    return if ( not $active_parameter_href->{vt_normalize} );

    my @to_process_references = check_references_for_vt(
        {
            active_parameter_href => $active_parameter_href,
            log                   => $log,
            parameter_href        => $parameter_href,
            vt_references_ref =>
              \@{ $active_parameter_href->{decompose_normalize_references} },
        }
    );

  REFERENCE:
    foreach my $reference_file_path (@to_process_references) {

        $log->info(q{[VT - Normalize and decompose]});
        $log->info( $TAB . q{File: } . $reference_file_path );

        ## Split multi allelic records into single records and normalize
        analysis_vt_core(
            {
                active_parameter_href   => $active_parameter_href,
                decompose               => 1,
                normalize               => 1,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                infile_path             => $reference_file_path,
                job_id_href             => $job_id_href,
                parameter_href          => $parameter_href,
            }
        );
    }
    return;
}

1;
