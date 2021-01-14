package MIP::Program::Gnu::Bash;

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
use MIP::Constants qw{ $NEWLINE $SINGLE_QUOTE $SPACE };
use MIP::List qw{ check_allowed_array_values };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {

    use base qw{ Exporter };
    require Exporter;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ gnu_cd gnu_export gnu_set gnu_trap gnu_ulimit gnu_unset gnu_wait };

}

sub gnu_cd {

##Function : Perl wrapper for writing cd recipe to already open $filehandle or return commands array. Based on cd 4.0
##Returns  : "@commands"
##Arguments: $directory_path         => Directory path
##         : $filehandle             => Filehandle to write to
##         : $stderrfile_path        => Stderrfile path
##         : $stderrfile_path_append => Append stderr info to file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $directory_path;
    my $stderrfile_path;
    my $stderrfile_path_append;

    my $tmpl = {
        directory_path => {
            store       => \$directory_path,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        stderrfile_path_append => {
            store       => \$stderrfile_path_append,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ cd };

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
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub gnu_export {

##Function : Perl wrapper for writing cd recipe to already open $filehandle or return commands array. Based on export 4.0
##Returns  : @commands
##Arguments: $bash_variable          => Bash variable to export
##         : $filehandle             => Filehandle to write to
##         : $stderrfile_path        => Stderrfile path
##         : $stderrfile_path_append => Append stderr info to file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $bash_variable;
    my $stderrfile_path;
    my $stderrfile_path_append;

    my $tmpl = {
        bash_variable => {
            required    => 1,
            store       => \$bash_variable,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        stderrfile_path_append => {
            store       => \$stderrfile_path_append,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ### cd
    ##Stores commands depending on input parameters
    my @commands = qw{ export };

    if ($bash_variable) {

        push @commands, $bash_variable;
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
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub gnu_set {

##Function : Perl wrapper for writing set recipe to already open $filehandle or return commands array. Based on set 4.0
##Returns  : @commands
##Arguments: $filehandle    => Filehandle to write to
##         : $separator     => Separator to use when writing
##         : $set_errexit   => Halt script if command has non-zero exit code (-e)
##         : $set_nounset   => Halt script if variable is uninitialised (-u)
##         : $set_pipefail  => Detect errors within pipes (-o pipefail)
##         : $unset_errexit => Unset errexit flag

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;

    ## Default(s)
    my $separator;
    my $set_errexit;
    my $set_nounset;
    my $set_pipefail;
    my $unset_errexit;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        set_errexit => {
            default     => 0,
            allow       => [ 0, 1 ],
            store       => \$set_errexit,
            strict_type => 1,
        },
        set_nounset => {
            default     => 0,
            allow       => [ 0, 1 ],
            store       => \$set_nounset,
            strict_type => 1,
        },
        set_pipefail => {
            default     => 0,
            allow       => [ 0, 1 ],
            store       => \$set_pipefail,
            strict_type => 1,
        },
        separator => {
            default     => $NEWLINE,
            store       => \$separator,
            strict_type => 1,
        },
        unset_errexit => {
            default     => 0,
            allow       => [ 0, 1 ],
            store       => \$unset_errexit,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ##Stores commands depending on input parameters
    my @commands;

    # Set flags
    if ($set_errexit) {

        push @commands, q{-e};
    }
    if ($set_nounset) {

        push @commands, q{-u};
    }
    if ($set_pipefail) {

        push @commands, q{-o pipefail};
    }
    if ($unset_errexit) {

        push @commands, q{+e};
    }

    ## Add set to each element
    @commands = map { q{set} . $SPACE . $_ } @commands;

    unix_write_to_file(
        {
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $NEWLINE,
        }
    );
    return @commands;
}

sub gnu_trap {

##Function : Perl wrapper for writing trap recipe to already open $filehandle or return commands array. Based on trap 4.0
##Returns  : @commands
##Arguments: $filehandle             => Filehandle to write to
##         : $stderrfile_path        => Stderrfile path
##         : $stderrfile_path_append => Append stderr info to file
##         : $trap_function_call     => Trap function argument
##         : $trap_signals_ref       => Array with signals to enable trap for {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $stderrfile_path;
    my $stderrfile_path_append;

    ## Default(s)
    my $trap_function_call;
    my $trap_signals_ref;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        stderrfile_path_append => {
            store       => \$stderrfile_path_append,
            strict_type => 1,
        },
        trap_function_call => {
            store       => \$trap_function_call,
            strict_type => 1,
        },
        trap_signals_ref => {
            allow => [
                sub {
                    check_allowed_array_values(
                        {
                            allowed_values_ref => [qw{ ERR EXIT TERM INT DEBUG }],
                            values_ref         => $arg_href->{trap_signals_ref},
                        }
                    );
                }
            ],
            default     => [],
            store       => \$trap_signals_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ trap };

    if ($trap_function_call) {

        # Quote function call to prevent word splitting
        push @commands, $SINGLE_QUOTE . $trap_function_call . $SINGLE_QUOTE;
    }
    if ( @{$trap_signals_ref} ) {

        push @commands, join $SPACE, @{$trap_signals_ref};
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
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub gnu_wait {

##Function : Perl wrapper for writing wait recipe to already open $filehandle or return commands array. Based on wait 4.0
##Returns  : @commands
##Arguments: $filehandle             => Filehandle to write to
##         : $processes_ref          => Specified processes to write to
##         : $stderrfile_path        => Stderrfile path
##         : $stderrfile_path_append => Append stderr info to file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $processes_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        processes_ref => {
            default     => [],
            store       => \$processes_ref,
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ wait };

    if ( @{$processes_ref} ) {

        push @commands, join $SPACE, @{$processes_ref};
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
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub gnu_ulimit {

## Function : Perl wrapper for setting bash ulimit to already open $filehandle or return commands array. Based on GNU Bash 4.0
## Returns  : @commands
## Arguments: $filehandle             => Filehandle to write to
##          : $max_open_files         => Set ulimit -n
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $max_open_files;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        max_open_files => {
            allow       => [ undef, qr/\A \d+ \z/xms ],
            store       => \$max_open_files,
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

    my @commands = qw{ ulimit };

    if ($max_open_files) {

        push @commands, q{-n} . $SPACE . $max_open_files;
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

sub gnu_unset {

## Function : Perl wrapper for writing unset recipe to already open $filehandle or return commands array. Based on unset 4.0
## Returns  : @commands
## Arguments: $bash_variable          => Variable to unset
##          : $filehandle             => Filehandle to write to
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bash_variable;
    my $filehandle;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        bash_variable => {
            defined     => 1,
            required    => 1,
            store       => \$bash_variable,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
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

    my @commands = qw{ unset };

    push @commands, $bash_variable;

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
