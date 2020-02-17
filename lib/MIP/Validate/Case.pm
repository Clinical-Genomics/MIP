package MIP::Validate::Case;

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
## MIPs lib/
use MIP::Constants qw{ $DOT $LOG_NAME $SINGLE_QUOTE $UNDERSCORE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ check_sample_ids };
}

sub check_sample_ids {

## Function : Check that the case_id and the sample_id(s) exists and are unique. Check if id sample_id contains "_".
## Returns  : 1
## Arguments: $case_id        => Case id
##          : $sample_ids_ref => Sample ids {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_id;
    my $sample_ids_ref;

    my $tmpl = {
      case_id => {
          defined     => 1,
          required    => 1,
          store       => \$case_id,
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

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Hash to test duplicate sample_ids later
    my %seen;

    if ( not @{$sample_ids_ref} ) {

        $log->fatal(q{Please provide sample_id(s)});
        exit 1;
    }

  SAMPLE_ID:
    foreach my $sample_id ( @{$sample_ids_ref} ) {

        ## Increment instance to check duplicates later
        $seen{$sample_id}++;

        ## Family_id cannot be the same as sample_id
        if ( $case_id eq $sample_id ) {

            $log->fatal( q{Case_id: '}
                  . $case_id
                  . q{' equals sample_id: '}
                  . $sample_id
                  . $SINGLE_QUOTE);
                  $log->fatal(q{Please make sure that the case_id and sample_id(s) are unique} );
            exit 1;
        }
        ## Check for unique sample_ids
        if ( $seen{$sample_id} > 1 ) {

            $log->fatal( q{Sample_id: '} . $sample_id . q{' is not uniqe.} );
            exit 1;
        }
        ## Sample_id contains "_", not allowed in filename convention
        if ( $sample_id =~ /$UNDERSCORE/sxm ) {

            $log->fatal( q{Sample_id: '}
                  . $sample_id
                  . q{' contains '_'});
                   $log->fatal(q{Please rename sample_id according to MIP's filename convention, removing the '_'.}
            );
            exit 1;
        }
    }
    return 1;
}

1;
