package MIP::Get::Executable;

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
use MIP::Constants qw{ $EMPTY_STR $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ get_executable };
}

sub get_executable {

## Function :
## Returns  :
## Arguments: $executable_name => Executable name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $executable_name;

    my $tmpl = {
        executable_name => {
            store       => \$executable_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %executable = ( vep =>
          { version_regexp => q?'if($_=~/ensembl-vep\s+:\s(\d+)/xms) {print $1;}'?, }, );

    if ( defined $executable_name and exists $executable{$executable_name} ) {

        return %{ $executable{$executable_name} };
    }
    return %executable;
}

1;
