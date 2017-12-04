package MIP::Check::Modules;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ check_perl_modules };
}

sub check_perl_modules {

## Function : Evaluate that all perl modules required by MIP are installed
## Returns  :
## Arguments: $modules_ref  => Array of module names
##          : $program_name => Program name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $modules_ref;
    my $program_name;

    my $tmpl = {
        modules_ref => {
            required    => 1,
            default     => [],
            strict_type => 1,
            store       => \$modules_ref
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    require Try::Tiny;
    use Try::Tiny;

  MODULE:
    foreach my $module ( @{$modules_ref} ) {

        ## Replace "::" with "/" since the automatic replacement magic only occurs for barewords.
        $module =~ s{::}{/}sxmg;

        ## Add perl module ending for the same reason
        $module .= q{.} . q{pm};

        try {
            require $module;
        }
        catch {
            croak(  q{NOTE: }
                  . $module
                  . q{ not installed - Please install to run }
                  . $program_name . qq{\n}
                  . q{NOTE: Aborting!} );
        };
    }
    return;
}

1;
