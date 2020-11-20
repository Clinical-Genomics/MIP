package MIP::Program::Conda;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use Getopt::Long;
use IO::Handle;
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use utf8;
use warnings qw{ FATAL utf8 };
use warnings;

## MIPs lib/
use MIP::Constants qw{ $LOG_NAME $NEWLINE $SPACE };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {

    use base qw{Exporter};
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.19;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ conda_activate conda_deactivate };
}

sub conda_activate {

##Function : Activate conda environment
##Returns  : @commands
##Arguments: $env_name   => Name of conda environment
##         : $filehandle => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $env_name;
    my $filehandle;

    my $tmpl = {
        env_name => {
            required    => 1,
            store       => \$env_name,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters

    # Basic command
    my @commands = qw{ conda activate };

    # Activates env, default base
    push @commands, $env_name;

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
        }
    );

    return @commands;
}

sub conda_deactivate {

##Function : Deactivate conda environment
##Returns  : @commands
##Arguments: $filehandle => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ conda deactivate };

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
        },
    );

    return @commands;
}

1;
