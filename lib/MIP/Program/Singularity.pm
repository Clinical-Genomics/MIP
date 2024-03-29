package MIP::Program::Singularity;

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
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ singularity_exec singularity_pull };
}

sub singularity_exec {

## Function : Perl wrapper for writing singularity execute command. Based on singularity v3.1.
## Returns  : @commands
## Arguments: $bind_paths_ref          => Array with paths to bind {REF}
##          : $container_cmds_ref      => Array with commands to be executed inside container {REF}
##          : $clean_env               => Start with clean environment
##          : $filehandle              => Filehandle to write to
##          : $gpu_switch              => Add nvidia experimental support
##          : $image                   => Singularity container name
##          : $no_home                 => Don't mount home if it isn't the current working dir
##          : $stderrfile_path         => Stderrfile path
##          : $stderrfile_path_append  => Append stderr info to file path
##          : $stdoutfile_path         => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bind_paths_ref;
    my $clean_env;
    my $container_cmds_ref;
    my $filehandle;
    my $gpu_switch;
    my $image;
    my $no_home;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        bind_paths_ref => {
            default     => [],
            store       => \$bind_paths_ref,
            strict_type => 1,
        },
        clean_env => {
            allow       => [ undef, 0, 1 ],
            store       => \$clean_env,
            strict_type => 1,
        },
        container_cmds_ref => {
            default     => [],
            store       => \$container_cmds_ref,
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
        gpu_switch => {
            allow       => [ undef, 0, 1 ],
            store       => \$gpu_switch,
            strict_type => 1,
        },
        no_home => {
            allow       => [ undef, 0, 1 ],
            store       => \$no_home,
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

    if ($gpu_switch) {
        push @commands, q{--nv};
    }
    ## Add bind paths
    if ( @{$bind_paths_ref} ) {
        push @commands, q{--bind} . $SPACE . join $COMMA, @{$bind_paths_ref};
    }

    if ($clean_env) {
        push @commands, q{--cleanenv};
    }

    if ($no_home) {
        push @commands, q{--no-home};
    }

    ## Add container
    push @commands, $image;

    ## Add optional commands to be run inside container
    if ( @{$container_cmds_ref} ) {
        push @commands, @{$container_cmds_ref};
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
            filehandle   => $filehandle,
            separator    => $SPACE,

        }
    );

    return @commands;
}

sub singularity_pull {

## Function : Perl wrapper for writing singularity pull command. Based on singularity v3.1.
## Returns  : @commands
## Arguments: $container_uri          => Container URI
##          : $filehandle             => Filehandle to write to
##          : $force                  => Force pull
##          : $outfile_path           => Save container to file
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $container_uri;
    my $filehandle;
    my $force;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        container_uri => {
            defined     => 1,
            required    => 1,
            store       => \$container_uri,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        force => {
            allow       => [ undef, 0, 1 ],
            store       => \$force,
            strict_type => 1,
        },
        outfile_path => {
            store       => \$outfile_path,
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
    my @commands = qw{ singularity pull };

    if ($force) {
        push @commands, q{--force};
    }

    if ($outfile_path) {
        push @commands, $outfile_path;
    }

    ## Add container uri
    push @commands, $container_uri;

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
1;
