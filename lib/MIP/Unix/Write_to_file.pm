package MIP::Unix::Write_to_file;

use strict;
use warnings;
use warnings qw(FATAL utf8);
use utf8;    #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );
use Carp;
use autodie;
use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;    #Do not convert to lower case

BEGIN {

    use base qw(Exporter);
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw(unix_write_to_file);
}

sub unix_write_to_file {

##unix_write_to_file

##Function : Perl wrapper for writing unix write to file recipe to already open $FILEHANDLE or return commands array.
##Returns  : ""
##Arguments: $commands_ref, $FILEHANDLE, $separator
##         : $commands_ref => Commands to write to file
##         : $FILEHANDLE => Filehandle to write to
##         : $separator  => Separator to use when writing

    my ($arg_href) = @_;

    ## Default(s)
    my $separator;

    ## Flatten argument(s)
    my $commands_ref;
    my $FILEHANDLE;

    my $tmpl = {
        commands_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$commands_ref
        },
        FILEHANDLE => { store => \$FILEHANDLE },
        separator  => {
            default     => q{ },
            strict_type => 1,
            store       => \$separator
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    ### Write to file helper module
    if ($FILEHANDLE) {

        print {$FILEHANDLE} join( $separator, @{$commands_ref} ) . $separator;
    }
    return;
}

1;
