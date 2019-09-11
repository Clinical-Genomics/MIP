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

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

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

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bind_paths_ref;
    my $FILEHANDLE;
    my $singularity_container;
    my $singularity_container_cmds_ref;

    ## Default(s)

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

    return @commands;
}

1;
