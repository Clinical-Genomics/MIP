package MIP::Language::Python;

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
use Readonly;

## MIPs lib/
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

## Constants
Readonly my $SPACE => q{ };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ python_core };

}

sub python_core {

## Function : Perl wrapper for writing python commands to $filehandle. Based on python 2.7/3.6.
## Returns  : @commands
## Arguments: $command_mode   => Execute python in command mode
##          : $module_mode    => Execute python in module mode
##          : $filehandle     => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $command_mode;
    my $filehandle;
    my $module_mode;

    my $tmpl = {
        command_mode => {
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$command_mode,
        },
        filehandle => {
            store => \$filehandle
        },
        module_mode => {
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$module_mode,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## command_mode and module_mode are mutually exclusive
    croak(q{Mutually exclusive options, $command_mode and $module_mode given})
      if ( $command_mode and $module_mode );

    # Stores commands depending on input parameters
    my @commands = q{python};

    ## Execute python command
    if ($command_mode) {
        push @commands, q{-c};
    }

    ## Execute python module
    if ($module_mode) {
        push @commands, q{-m};
    }

    unix_write_to_file(
        {
            commands_ref => \@commands,
            separator    => $SPACE,
            filehandle   => $filehandle,
        }
    );
    return @commands;
}

1;
