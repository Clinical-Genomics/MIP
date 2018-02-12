package MIP::Package_manager::Cpanm;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;

## CPANM
use Readonly;

## MIPs lib/
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

## Constants
Readonly my $SPACE => q{ };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.0.1;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ cpanm_install };

}

sub cpanm_install {

## Function  : Perl wrapper for writing cpanm recipe to $FILEHANDLE.
## Returns   :
## Arguments : $FILEHANDLE  => Filehandle to write to
##           : $force       => Force install
##           : $modules_ref => Perl modules {REF}
##           : $quiet       => Supress output

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $force;
    my $modules_ref;
    my $quiet;

    my $tmpl = {
        modules_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$modules_ref,
            strict_type => 1,
        },
        FILEHANDLE => {
            required => 1,
            store    => \$FILEHANDLE,
        },
        force => {
            allow       => [ 0, 1 ],
            default     => 1,
            store       => \$force,
            strict_type => 1,
        },
        quiet => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$quiet,
            strict_type => 1,
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
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );
    return @commands;
}

1;
