package MIP::Check::Vcfparser;

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
use MIP::Constants qw{ $LOG $NEWLINE $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ check_vcfparser_cli };
}

sub check_vcfparser_cli {

## Function : Check basic vcfparser CLI options
## Returns  :
## Arguments: $range_feature_annotation_column => Range feature file annotation columns
##          : $range_feature_file              => Range feature file
##          : $select_feature_file             => Select feature file
##          : $select_outfile             => Select file path
##          : $select_feature_matching_column  => Select feature matching column

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $range_feature_annotation_column;
    my $range_feature_file;
    my $select_outfile;
    my $select_feature_matching_column;

    ## Default(s)
    my $select_feature_file;

    my $tmpl = {
        range_feature_annotation_column => {
            store       => \$range_feature_annotation_column,
            strict_type => 1,
        },
        range_feature_file => {
            store       => \$range_feature_file,
            strict_type => 1,
        },
        select_feature_file => {
            default     => 0,
            store       => \$select_feature_file,
            strict_type => 1,
        },
        select_feature_matching_column => {
            store       => \$select_feature_matching_column,
            strict_type => 1,
        },
        select_outfile => { store => \$select_outfile, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Log::MIP_log4perl qw{ retrieve_log };

    my $log = retrieve_log( { log_name => $LOG, } );

    my %cli = (
        range_feature_annotation_column => {
            dependens_on_value => $range_feature_file,
            msg =>
              q{Need to specify which feature column(s) to use with range feature file: }
              . $range_feature_file
              . q{ when annotating variants by using flag -rf_ac}
              . $NEWLINE,
            value => $range_feature_annotation_column,
        },
        select_feature_matching_column => {
            dependens_on_value => $select_feature_file,
            msg =>
              q{Need to specify which feature column to use with select feature file: }
              . $select_feature_file
              . q{ when selecting variants by using flag -sf_mc}
              . $NEWLINE,
            value => $select_feature_matching_column,
        },
        select_outfile => {
            dependens_on_value => $select_feature_file,
            msg =>
q{Need to specify which a select outfile to use when selecting variants by using flag -sof}
              . $NEWLINE,
            value => $select_outfile,
        }
    );

  OPTION:
    foreach my $option ( keys %cli ) {

        my $option_value      = $cli{$option}{value};
        my $option_dependency = $cli{$option}{dependens_on_value};

        if ( not $option_value and $option_dependency ) {

            $log->fatal( $cli{$option}{msg} );
            exit 1;
        }
    }
    return;
}

1;
