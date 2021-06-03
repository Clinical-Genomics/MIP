package MIP::Unix::Write_to_file;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

## MIPs lib/
use MIP::Constants qw{ $NEWLINE $SPACE };

BEGIN {

    use base qw{ Exporter };
    require Exporter;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ unix_write_to_file };
}

sub unix_write_to_file {

## Function : Perl wrapper for writing unix write to file recipe to already open $filehandle.
## Returns  :
## Arguments: $commands_ref => Commands to write to file
##          : $filehandle   => Filehandle to write to
##          : $separator    => Separator to use when writing

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $commands_ref;
    my $filehandle;

    ## Default(s)
    my $separator;

    my $tmpl = {
        commands_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$commands_ref,
            strict_type => 1,
        },
        filehandle => { store => \$filehandle, },
        separator  => {
            default     => q{ },
            store       => \$separator,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ### Write to file helper module
    if ( $filehandle && @{$commands_ref} ) {

        ## Write each element in line
        if ( $separator ne q{\n} ) {

            print {$filehandle} join( $separator, @{$commands_ref} ) . $separator;
        }
        else {
            ## Write each command per line

          LINES:
            foreach my $line ( @{$commands_ref} ) {

                say {$filehandle} $line;
            }
            print {$filehandle} $NEWLINE;
        }
    }
    return;
}

1;
