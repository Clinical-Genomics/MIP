package MIP::Check::Parameter;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie;
use Email::Valid;
use Readonly;
use List::MoreUtils qw { any uniq };

## MIPs lib/
use MIP::Constants qw{ $LOG_NAME };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.48;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      check_active_installation_parameters
    };
}

sub check_active_installation_parameters {

## Function : Some active_parameter checks that are common to both installations. Returns "1" if all is OK
## Returns  : 1 or exit
## Arguments: $project_id => Project id
##          : sbatch_mode => Sbatch mode boolean

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $project_id;
    my $sbatch_mode;

    my $tmpl = {
        project_id => {
            store       => \$project_id,
            strict_type => 1,
        },
        sbatch_mode => {
            store       => \$sbatch_mode,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Check that a project id has been set if SBATCH mode
    if ( $sbatch_mode
        and not $project_id )
    {
        $log->fatal(
q{The parameter "project_id" must be set when a sbatch installation has been requested}
        );
        exit 1;
    }
    return 1;
}

1;
