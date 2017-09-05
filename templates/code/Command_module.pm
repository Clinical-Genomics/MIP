#package MIP::PATH::TO::MODULE                                # adjust code here

use strict;
use warnings;
use warnings qw( FATAL utf8 );
use utf8;    #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );
use Carp;
use English qw{-no_match_vars};
use Params::Check qw[check allow last_error];

use Readonly;

use FindBin qw{$Bin};    #Find directory of script
use File::Basename qw{dirname};
use File::Spec::Functions qw{catdir};

## MIPs lib/
use lib catdir( dirname($Bin), 'lib' );
use MIP::Unix::Standard_streams qw{unix_standard_streams};
use MIP::Unix::Write_to_file qw{unix_write_to_file};

BEGIN {
    require Exporter;
    use base qw(Exporter);

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw(MIP::PATH::TO::MODULE);               # adjust code here
}

## Constants (add as many as you need)
Readonly my $SPACE => q{ };


sub name_of_subroutine {                                      # adjust code here

##bwa_mem

##Function : Perl wrapper for generic commands module.
##Returns  : "@commands"
##Arguments: $stdoutfile_path, $stderrfile_path, stderrfile_path_append, $FILEHANDLE
##         : $stdoutfile_path        => Stdoutfile path
##         : $stderrfile_path        => Stderrfile path
##         : $stderrfile_path_append => Append stderr info to file path
##         : $FILEHANDLE             => Filehandle to write to

    my ($arg_href) = @_;

    ## Default(s)

    ## Flatten argument(s)
    my $stdoutfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $FILEHANDLE;

    # initialize other command arguments here
    my $array_argument;
    my $scalar_argument;

    my $tmpl = {
        stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
        FILEHANDLE    => { store => \$FILEHANDLE },

        # add all the above command arguments here, following this specification:
        #array_argument => {
        #  default     => qw{ option1 },
        #  allow       => [qw{option1 option2 option3}],
        #  strict_type => 1,
        #  store       => \$array_argument
        #},

        #scalar_argument => { strict_type => 1, store => \$scalar_argument },
    };

    check( $tmpl, $arg_href, 1 ) or croak qw{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw(BASE_COMMAND);

    ##Options (add one for each argument specified above).
    #if ($argument_provided) {

        # Add argument to commands
        #push @commands, q{>} . $SPACE . $argument_provided;
    #}

    push @commands,
      unix_standard_streams(
        {
            stdoutfile_path        => $stdoutfile_path,
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
