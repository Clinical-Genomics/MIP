package MIP::Star;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

## MIPs lib/
use MIP::Constants qw{ $LOG_NAME };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ check_interleaved_files_for_star };
}

sub check_interleaved_files_for_star {

## Function : Check for interleaved fastq infiles. Log and exit if true.
## Returns  :
## Arguments: $file_info_href          => File info hash {REF}
##          : $sample_ids_ref          => Sample ids {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $sample_ids_ref;

    my $tmpl = {
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        sample_ids_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$sample_ids_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File_info qw{ get_sample_file_attribute };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

  SAMPLE_ID:
    foreach my $sample_id ( @{$sample_ids_ref} ) {

        my %attribute = get_sample_file_attribute(
            {
                file_info_href => $file_info_href,
                sample_id      => $sample_id,
            }
        );

      INFILES:
        foreach my $file_name ( @{ $attribute{mip_infiles} } ) {

            my $is_interleaved = get_sample_file_attribute(
                {
                    attribute      => q{is_interleaved},
                    file_info_href => $file_info_href,
                    file_name      => $file_name,
                    sample_id      => $sample_id,
                }
            );

            ## STAR does not support interleaved fastq files
            if ($is_interleaved) {
                $log->fatal(q{MIP rd_rna does not support interleaved fastq files});
                $log->fatal( q{Please deinterleave: }
                      . catfile( $attribute{mip_infiles_dir}, $file_name ) );
                exit 1;
            }
        }
    }
    return;
}

1;
