package MIP::Program::Docker;

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
use MIP::Constants qw{ $COMMA $DOUBLE_QUOTE $SPACE };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ docker_pull docker_run };
}

sub docker_pull {

## Function : Perl wrapper for pulling docker images Based on Docker 19.03.8
## Returns  : @commands
## Arguments: $filehandle             => Filehandle to write to
##          : $image                  => Image to pull
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $image;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        image => {
            defined     => 1,
            required    => 1,
            store       => \$image,
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
    my @commands = qw{ docker pull };

    push @commands, $image;

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
            filehandle   => $filehandle,
            separator    => $SPACE,

        }
    );
    return @commands;
}

sub docker_run {

## Function : Perl wrapper for docker run command. Based on Docker 19.03.8
## Returns  : @commands
## Arguments: $bind_paths_ref         => Bind host directory to container {REF}
##          : $container_cmds_ref     => Cmds to be executed in container {REF}
##          : $entrypoint             => Override container entrypoint
##          : $filehandle             => Filehandle to write to
##          : $image                  => Image to run
##          : $remove                 => Remove stopped container
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdinfile_path         => Stdinfile path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bind_paths_ref;
    my $container_cmds_ref;
    my $filehandle;
    my $image;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdinfile_path;
    my $stdoutfile_path;

    ## Default(s)
    my $entrypoint;
    my $remove;

    my $tmpl = {
        bind_paths_ref => {
            default     => [],
            store       => \$bind_paths_ref,
            strict_type => 1,
        },
        container_cmds_ref => {
            default     => [],
            store       => \$container_cmds_ref,
            strict_type => 1,
        },
        entrypoint => {
            default     => $DOUBLE_QUOTE . $DOUBLE_QUOTE,
            store       => \$entrypoint,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        image => {
            defined     => 1,
            required    => 1,
            store       => \$image,
            strict_type => 1,
        },
        remove => {
            allow       => [ 0, 1 ],
            default     => 1,
            store       => \$remove,
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
        stdinfile_path  => { store => \$stdinfile_path, strict_type => 1, },
        stdoutfile_path => {
            store       => \$stdoutfile_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ docker run };

    if ( @{$bind_paths_ref} ) {

        push @commands, map { q{--volume} . $SPACE . $_ } @{$bind_paths_ref};
    }

    if ($remove) {

        push @commands, q{--rm};
    }

    if ($entrypoint) {

        push @commands, q{--entrypoint} . $SPACE . $entrypoint;
    }

    push @commands, $image;

    if ($container_cmds_ref) {

        push @commands, @{$container_cmds_ref};
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
            stdinfile_path         => $stdinfile_path,
            stdoutfile_path        => $stdoutfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,

        }
    );
    return @commands;
}

1;
