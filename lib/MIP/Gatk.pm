package MIP::Gatk;

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

## MIPs lib
use MIP::Constants qw{ $LOG_NAME $NEWLINE $TAB };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ check_gatk_sample_map_paths };
}

sub check_gatk_sample_map_paths {

## Function : Check that the supplied gatk sample map file paths exists
## Returns  :
## Arguments: $sample_map_path => Sample map path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_map_path;

    my $tmpl = {
        sample_map_path => {
            defined     => 1,
            required    => 1,
            store       => \$sample_map_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Io::Read qw{ read_from_file };

    ## Constants
    Readonly my $FIELD_COUNTER => 2;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my @missing_files;

    my @lines = read_from_file(
        {
            chomp  => 1,
            format => q{line_by_line},
            path   => $sample_map_path,
        }
    );

  LINE:
    foreach my $line (@lines) {

        ## Get sample and file path (and check proper format)
        my ( $sample, $file_path, $unexpected_data ) =
          split $TAB, $line, $FIELD_COUNTER + 1;

        ## Make sure that we get what we expect
        if ( defined $unexpected_data ) {

            $log->logcarp(
                q{Unexpected trailing garbage at end of line '} . $line . q{':},
                $NEWLINE . $TAB . $unexpected_data . $NEWLINE );
        }

        ## Path exists
        next LINE if ( -e $file_path );

        my $error_msg =
            q{The supplied file path: }
          . $file_path
          . q{ from sample map file: }
          . $sample_map_path
          . q{ does not exist};
        push @missing_files, $error_msg;
    }
    return 1 if ( not @missing_files );

    ## Broadcast missing files
    foreach my $error_msg (@missing_files) {
        $log->fatal($error_msg);
    }
    exit 1;
}

1;
