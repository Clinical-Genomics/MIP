package MIP::Get::Analysis;

use strict;
use warnings;
use warnings qw( FATAL utf8 );
use utf8;    # Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );
use Carp;
use autodie;
use Params::Check qw[check allow last_error];

use FindBin qw($Bin);    # Find directory of script
use File::Basename qw(dirname);
use File::Spec::Functions qw(catdir);
use Readonly;

## Third party module(s)
use List::Util qw(all);

## MIPs lib/
use lib catdir( dirname($Bin), 'lib' );

BEGIN {
    use base qw (Exporter);
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw(get_overall_analysis_type);

}

sub get_overall_analysis_type {

##get_overall_analysis_type

##Function : Detect if all samples has the same sequencing type and return consensus or mixed
##Returns  : "consensus/mixed analysis_type"
##Arguments: $analysis_type_hef
##         : $analysis_type_hef => The analysis_type hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $analysis_type_hef;

    my $tmpl = {
        analysis_type_hef => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$analysis_type_hef
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw[Could not parse arguments!];

    my @analysis_types = (qw{wes wgs rapid});

    foreach my $analysis_type (@analysis_types) {

        if ( all { $_ eq $analysis_type } values %{$analysis_type_hef} ) {

            return $analysis_type;
        }
    }

    # No consensus, then it must be mixed
    return q{mixed};
}

1;
