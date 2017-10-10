package MIP::Check::Modules_v5_10;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use English qw{ -no_match_vars };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ check_perl_modules };
}

sub check_perl_modules {

## Function : Evaluate whether the perl modules are installed and returns any missing modules.
## Returns  : @missing_modules
##          : $modules_ref  => Array of module names {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $modules_ref = $arg_href->{modules_ref};

    my $module_file;
    my @missing_modules;

  MODULE:
    foreach my $module ( @{$modules_ref} ) {

        ## Add perl module ending since the automatic replacement magic only occurs for barewords.
        $module_file = $module . q{.} . q{pm};

        ## Replace "::" with "/" for the same reason.
        $module_file =~ s/::/\//sxmg;

        eval { require $module_file; };
        if ($EVAL_ERROR) {
            push @missing_modules, $module;
        }
    }
    return @missing_modules;
}

1;
