package MIP::Recipes::Build::Cancer;

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

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ build_cancer_meta_files };
}

## Constants
Readonly my $SPACE     => q{ };
Readonly my $EMPTY_STR => q{};
Readonly my $TAB       => qq{\t};

sub build_cancer_meta_files {

## Function : Pipeline recipe for cancer data analysis.
## Returns  :

## Arguments: $parameter_href          => Parameter hash {REF}
##          : $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $log                     => Log object to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $log;

    my $tmpl = {
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Reference
      qw{  check_parameter_metafiles check_references_for_vt };

    $log->info( $TAB . q{Reference check: Reference prerequisites checked} );

## Check if vt has processed references, if not try to reprocesses them before launcing modules
    $log->info(q{[Reference check - Reference processed by VT]});
    if (   $active_parameter_href->{vt_decompose}
        || $active_parameter_href->{vt_normalize} )
    {

        my @to_process_references = check_references_for_vt(
            {
                parameter_href        => $parameter_href,
                active_parameter_href => $active_parameter_href,
                vt_references_ref =>
                  \@{ $active_parameter_href->{decompose_normalize_references}
                  },
                log => $log,
            }
        );

      REFERENCE:
        foreach my $reference_file_path (@to_process_references) {

            $log->info(q{[VT - Normalize and decompose]});
            $log->info( $TAB . q{File: } . $reference_file_path );

            ## Split multi allelic records into single records and normalize
            analysis_vt_core(
                {
                    parameter_href          => $parameter_href,
                    active_parameter_href   => $active_parameter_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    infile_path             => $reference_file_path,
                    program_directory       => q{vt},
                    decompose               => 1,
                    normalize               => 1,
                }
            );
        }
    }
    return;
}

1;
