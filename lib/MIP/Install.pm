package MIP::Install;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
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
    our @EXPORT_OK = qw{ check_mip_executable };
}

sub check_mip_executable {

## Function : Check if mip installation exists and is executable
## Returns  :
##          : $conda_prefix_path => Conda prefix path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_prefix_path;

    my $tmpl = {
        conda_prefix_path => {
            defined     => 1,
            required    => 1,
            store       => \$conda_prefix_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $log = Log::Log4perl->get_logger($LOG_NAME);

    return 1 if ( not -x catfile( $conda_prefix_path, qw{ bin mip } ) );

    $log->info(q{MIP is already installed in the specified conda environment.});

    $log->warn(q{This will overwrite the current installation of MIP});
    return;
}

1;
