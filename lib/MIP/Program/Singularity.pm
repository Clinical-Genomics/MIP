package MIP::Program::Singularity;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ singularity_exec };
}

sub singularity_exec {

## Function : Perl wrapper for writing singularity execute command. Based on singularity v3.1.
## Returns  : @commands
## Arguments: $bind_paths_ref                 => Array with paths to bind {REF}
##          : $FILEHANDLE                     => Filehandle to write to
##          : $singularity_container          => Singularity container name
##          : $singularity_container_cmds_ref => Array with commands to be executed inside container {REF} 
##          : $stderrfile_path                => Stderrfile path
##          : $stderrfile_path_append         => Append stderr info to file path
##          : $stdoutfile_path                => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bind_paths_ref;
    my $FILEHANDLE;
    my $singularity_container;
    my $singularity_container_cmds_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        bind_paths_ref => {
            default     => [],
            store       => \$bind_paths_ref,
            strict_type => 1,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        singularity_container => {
            defined     => 1,
            required    => 1,
            store       => \$singularity_container,
            strict_type => 1,
        },
        singularity_container_cmds_ref => {
            default     => [],
            store       => \$singularity_container_cmds_ref,
            strict_type => 1,
        },
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        stderrfile_path_append => {
            store       => \$stderrfile_path_append,
            strict_type => 1,
        },
        stdoutfile_path => {
            store       => \$stdoutfile_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ singularity exec };

    ## Add bind paths
    if ( @{$bind_paths_ref} ) {
        push @commands, q{--bind} . $SPACE . join $COMMA, @{$bind_paths_ref};
    }

    ## Add container
    push @commands, $singularity_container;

    ## Add optional commands to be run inside container
    if ( @{$singularity_container_cmds_ref} ) {
        push @commands, @{$singularity_container_cmds_ref};
    }
    
    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
            stdoutfile_path        => $stdoutfile_path,
        }
      );

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
