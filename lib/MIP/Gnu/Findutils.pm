package MIP::Gnu::Findutils;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;    #Allow unicode characters in this script
use open qw{  :encoding(UTF-8) :std};
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };

use Readonly;

use FindBin qw{ $Bin };    #Find directory of script
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir };

## MIPs lib/
use lib catdir( dirname($Bin), q{ lib } );
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ xargs };
}

## Constants
Readonly my $DOUBLE_QUOTE => q{"};
Readonly my $SPACE        => q{ };

sub xargs {

##xargs

##Function : Perl wrapper for writing xargs recipe to already open $FILEHANDLE or return command line string. Based on xargs 4.4.2
##Returns  : "@commands"
##Arguments: $shell_commands_ref, $FILEHANDLE, $replace_str, $verbose, $max_args, $max_procs, $placeholder_symbol
##         : $shell_commands_ref => The string following this command will be interpreted as a shell command {REF}
##         : $FILEHANDLE         => Filehandle to write to
##         : $replace_str        => Replace string.  Enables us to tell xargs where to put the command file lines
##         : $verbose            => Print the command line on the standard error output before executing it
##         : $max_args           => Use at most max-args arguments per command line
##         : $max_procs          => Run up to max-procs processes at a time
##         : $placeholder_symbol => Set placeholder symbol

    my ($arg_href) = @_;

    ## Default(s)
    my $replace_str;
    my $verbose;
    my $max_args;
    my $max_procs;
    my $placeholder_symbol;

    ## Flatten argument(s)
    my $shell_commands_ref;
    my $FILEHANDLE;

    my $tmpl = {
        shell_commands_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$shell_commands_ref
        },
        FILEHANDLE  => { store => \$FILEHANDLE },
        replace_str => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$replace_str
        },
        verbose => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$verbose
        },
        max_args => {
            default     => 1,
            allow       => qr/ ^\d+$ /sxm,
            strict_type => 1,
            store       => \$max_args
        },
        max_procs => {
            default     => 1,
            allow       => qr/ ^\d+$ /sxm,
            strict_type => 1,
            store       => \$max_procs
        },
        placeholder_symbol => {
            default     => q?{}?,
            strict_type => 1,
            store       => \$placeholder_symbol
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Xargs
    # Stores commands depending on input parameters
    my @commands = qw( xargs );

    if ($replace_str) {

     # Replace-str; Enables us to tell xargs where to put the command file lines
        push @commands, q{-i};
    }
    if ($verbose) {

       # Print the command line on the standard error output before executing it
        push @commands, q{--verbose};
    }
    if ($max_args) {

        # Use at most max-args arguments per command line
        push @commands, q{-n} . $SPACE . $max_args;
    }
    if ($max_procs) {

        # Run up to max-procs processes at a time
        push @commands, q{-P} . $SPACE . $max_procs;
    }
    if ($placeholder_symbol) {

      # The string following this command will be interpreted as a shell command
        push @commands, q{sh -c} . $SPACE . $DOUBLE_QUOTE;

        if ($shell_commands_ref) {

            push @commands, join $SPACE, @{$shell_commands_ref};
        }

        # Set placeholder and end quotes
        push @commands, $placeholder_symbol . $DOUBLE_QUOTE . $SPACE;
    }
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
