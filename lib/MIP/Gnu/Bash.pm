package MIP::Gnu::Bash;

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

use FindBin qw($Bin);                 #Find directory of script
use File::Basename qw(dirname);
use File::Spec::Functions qw(catdir);

## MIPs lib/
use lib catdir( dirname($Bin), 'lib' );
use MIP::Unix::Standard_streams qw(unix_standard_streams);
use MIP::Unix::Write_to_file qw(unix_write_to_file);

BEGIN {

    use base qw(Exporter);
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw(gnu_cd gnu_trap);

}

## Constants
my $SPACE      = q{ };
my $APOSTROPHE = q{'};

sub gnu_cd {

##gnu_cd

##Function : Perl wrapper for writing cd recipe to already open $FILEHANDLE or return commands array. Based on cd 4.0
##Returns  : "@commands"
##Arguments: $FILEHANDLE, $directory_path, $stderrfile_path, $stderrfile_path_append
##         : $FILEHANDLE             => Filehandle to write to
##         : $directory_path         => Directory path
##         : $stderrfile_path        => Stderrfile path
##         : $stderrfile_path_append => Append stderr info to file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $directory_path;
    my $stderrfile_path;
    my $stderrfile_path_append;

    my $tmpl = {
        FILEHANDLE      => { store       => \$FILEHANDLE },
        directory_path  => { strict_type => 1, store => \$directory_path },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    ### cd
    ##Stores commands depending on input parameters
    my @commands = qw(cd);

    ## Options
    if ($directory_path) {

        push @commands, $directory_path;
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
        }
      );
    unix_write_to_file(
        {
            commands_ref => \@commands,
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    return @commands;
}

sub gnu_trap {

##gnu_trap

##Function : Perl wrapper for writing trap recipe to already open $FILEHANDLE or return commands array. Based on trap 4.0
##Returns  : "@commands"
##Arguments: $trap_signals_ref, $FILEHANDLE, $stderrfile_path, $stderrfile_path_append, $trap_function_call
##         : $trap_signals_ref       => Array with signals to enable trap for {REF}
##         : $trap_function_call     => The trap function argument
##         : $FILEHANDLE             => Filehandle to write to
##         : $stderrfile_path        => Stderrfile path
##         : $stderrfile_path_append => Append stderr info to file

    my ($arg_href) = @_;

    ## Default(s)
    my $trap_signals_ref;
    my $trap_function_call;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $stderrfile_path;
    my $stderrfile_path_append;

    use MIP::Check::Parameter qw(check_allowed_array_values);

    my $tmpl = {
        trap_signals_ref => {
            default => [],
            allow   => [
                sub {
                    check_allowed_array_values(
                        {
                            allowed_values_ref => [qw(ERR EXIT TERM INT DEBUG)],
                            values_ref         => $arg_href->{trap_signals_ref},
                        }
                    );
                }
            ],
            strict_type => 1,
            store       => \$trap_signals_ref
        },
        FILEHANDLE      => { store       => \$FILEHANDLE },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
        trap_function_call => {
            strict_type => 1,
            store       => \$trap_function_call
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    ### trap
    ##Stores commands depending on input parameters
    my @commands = qw(trap);

    ## Options

    if ($trap_function_call) {

        # Quote function call to prevent word splitting
        push @commands, $APOSTROPHE . $trap_function_call . $APOSTROPHE;
    }
    if ( @{$trap_signals_ref} ) {

        push @commands, join( $SPACE, @{$trap_signals_ref} );
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
        }
      );
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
