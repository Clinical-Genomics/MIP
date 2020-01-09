package MIP::File::Path;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd qw{ abs_path };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ get_absolute_path };
}

sub get_absolute_path {

## Function : Get absolute path for supplied path or croaks and exists if path does not exists
## Returns  : $path (absolute path)
## Arguments: $parameter_name => Parameter to be evaluated
##          : $path           => Supplied path to be updated/evaluated

    my ($arg_href) = @_;

    ##Flatten argument(s)
    my $parameter_name;
    my $path;

    my $tmpl = {
        parameter_name => {
            defined     => 1,
            required    => 1,
            store       => \$parameter_name,
            strict_type => 1,
        },
        path => { defined => 1, required => 1, store => \$path, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Reformat to absolute path
    my $absolute_path = abs_path($path);

    return $absolute_path if ( defined $absolute_path );

    croak(  q{Could not find absolute path for }
          . $parameter_name . q{: }
          . $path
          . q{. Please check the supplied path!} );
}

1;
