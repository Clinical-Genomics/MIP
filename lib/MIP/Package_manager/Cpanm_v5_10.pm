package MIP::Package_manager::Cpanm_v5_10;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;    #Allow unicode characters in this script
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.0.1;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ cpanm_install_module };

}

sub cpanm_install_module {

## Function  : Perl wrapper for writing cpanm recipe to array (@commands).
## Returns   : @commands
##           : $force       => Force install
##           : $modules_ref => Perl modules {REF}
##           : $notest      => No testing of installed cpanmodules
##           : $quiet       => Supress output
##           : $verbose     => Toggle verbsoe output

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $force       = $arg_href->{force};
    my $modules_ref = $arg_href->{modules_ref};
    my $notest      = $arg_href->{notest};
    my $quiet       = $arg_href->{quiet};
    my $verbose     = $arg_href->{verbose};

    ## Base command
    my @commands = q{cpanm};

    ## Add optional force flag
    if ($force) {
        push @commands, q{--force};
    }

    ## Add optional quiet flag
    if ($quiet) {
        push @commands, q{--quiet};
    }

    ## Add optional verbose flag
    if ($verbose) {
        push @commands, q{--verbose};
    }

    ## Add optional no test flag
    if ($notest) {
        push @commands, q{--notest};
    }

    ## Add cpan modules
    push @commands, join q{ }, @{$modules_ref};

    return @commands;
}

1;
