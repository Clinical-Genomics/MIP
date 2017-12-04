package MIP::Package_manager::Cpanm;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;    #Allow unicode characters in this script
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use Readonly;

## MIPs lib/
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

## Constants
Readonly my $SPACE => q{ };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.0.0;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ cpanm_install };

}

sub cpanm_install {

## cpanm_install

## Function  : Perl wrapper for writing cpanm recipe to $FILEHANDLE.
## Returns   : ""
## Arguments : $modules_ref, $FILEHANDLE, $force, $quiet
##           : $modules_ref => Perl modules {REF}
##           : $FILEHANDLE  => Filehandle to write to
##           : $force       => Force install
##           : $quiet       => Supress output

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $modules_ref;
    my $FILEHANDLE;
    my $force;
    my $quiet;

    my $tmpl = {
        modules_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$modules_ref
        },
        FILEHANDLE => {
            required => 1,
            store    => \$FILEHANDLE
        },
        force => {
            default     => 1,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$force,
        },
        quiet => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$quiet,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

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

    ## Add cpan modules
    push @commands, join $SPACE, @{$modules_ref};

    unix_write_to_file(
        {
            commands_ref => \@commands,
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;
}

1;
