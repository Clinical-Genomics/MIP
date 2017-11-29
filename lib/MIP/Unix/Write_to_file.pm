package MIP::Unix::Write_to_file;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie;
use Readonly;

BEGIN {

    use base qw(Exporter);
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ unix_write_to_file };
}

## Constants
Readonly my $NEWLINE => qq{\n};

sub unix_write_to_file {

## Function : Perl wrapper for writing unix write to file recipe to already open $FILEHANDLE.
## Returns  :
## Arguments: $commands_ref => Commands to write to file
##          : $FILEHANDLE   => Filehandle to write to
##          : $separator    => Separator to use when writing

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $commands_ref;
    my $FILEHANDLE;

    ## Default(s)
    my $separator;

    my $tmpl = {
        commands_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$commands_ref,
        },
        FILEHANDLE => { store => \$FILEHANDLE, },
        separator  => {
            default     => q{ },
            strict_type => 1,
            store       => \$separator,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ### Write to file helper module
    if ( $FILEHANDLE && @{$commands_ref} ) {

        ## Write each element in line
        if ( $separator ne q{\n} ) {

            print {$FILEHANDLE} join( $separator, @{$commands_ref} )
              . $separator;
        }
        else {
            ## Write each command per line

          LINES:
            foreach my $line ( @{$commands_ref} ) {

                say {$FILEHANDLE} $line;
            }
            print {$FILEHANDLE} $NEWLINE;
        }
    }
    return;
}

1;
