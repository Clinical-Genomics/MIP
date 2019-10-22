package MIP::Gnu::Bash;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use autodie;
use Params::Check qw{ check allow last_error };

## MIPs lib/
use MIP::Constants qw{ $NEWLINE $SINGLE_QUOTE $SPACE };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {

    use base qw{ Exporter };
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ gnu_cd gnu_trap gnu_set gnu_unset gnu_wait };

}

sub gnu_cd {

##Function : Perl wrapper for writing cd recipe to already open $FILEHANDLE or return commands array. Based on cd 4.0
##Returns  : "@commands"
##Arguments: $directory_path         => Directory path
##         : $FILEHANDLE             => Filehandle to write to
##         : $stderrfile_path        => Stderrfile path
##         : $stderrfile_path_append => Append stderr info to file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $directory_path;
    my $stderrfile_path;
    my $stderrfile_path_append;

    my $tmpl = {
        directory_path => {
            strict_type => 1,
            store       => \$directory_path,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path,
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ### cd
    ##Stores commands depending on input parameters
    my @commands = q{cd};

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
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub gnu_trap {

##Function : Perl wrapper for writing trap recipe to already open $FILEHANDLE or return commands array. Based on trap 4.0
##Returns  : "@commands"
##Arguments: $FILEHANDLE             => Filehandle to write to
##         : $stderrfile_path        => Stderrfile path
##         : $stderrfile_path_append => Append stderr info to file
##         : $trap_function_call     => The trap function argument
##         : $trap_signals_ref       => Array with signals to enable trap for {REF}

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
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path,
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append,
        },
        trap_function_call => {
            strict_type => 1,
            store       => \$trap_function_call,
        },
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
            store       => \$trap_signals_ref,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ### trap
    ##Stores commands depending on input parameters
    my @commands = q{trap};

    ## Options

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
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    return @commands;
}

sub gnu_set {

##Function : Perl wrapper for writing set recipe to already open $FILEHANDLE or return commands array. Based on set 4.0
##Returns  : "@commands"
##Arguments: $FILEHANDLE             => Filehandle to write to
##         : $separator              => Separator to use when writing
##         : $set_errexit            => Halt script if command has non-zero exit code (-e)
##         : $set_nounset            => Halt script if variable is uninitialised (-u)
##         : $set_pipefail           => Detect errors within pipes (-o pipefail)
##         : $unset_errexit          => Unset errexit flag

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;

    ## Default(s)
    my $separator;
    my $set_errexit;
    my $set_nounset;
    my $set_pipefail;
    my $unset_errexit;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
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
            FILEHANDLE   => $FILEHANDLE,
            separator    => $NEWLINE,
        }
    );
    return @commands;
}

sub gnu_wait {

##Function : Perl wrapper for writing wait recipe to already open $FILEHANDLE or return commands array. Based on wait 4.0
##Returns  : "@commands"
##Arguments: $FILEHANDLE             => Filehandle to write to
##         : $processes_ref          => Specified processes to write to
##         : $stderrfile_path        => Stderrfile path
##         : $stderrfile_path_append => Append stderr info to file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $processes_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        processes_ref => {
            default     => [],
            strict_type => 1,
            store       => \$processes_ref,
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path,
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ### wait
    ##Stores commands depending on input parameters
    my @commands = q{wait};

    ## Options
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
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    return @commands;
}

sub gnu_unset {

## Function : Perl wrapper for writing unset recipe to already open $FILEHANDLE or return commands array. Based on unset 4.0
## Returns  : @commands
## Arguments: $bash_variable          => Variable to unset
##          : $FILEHANDLE             => Filehandle to write to
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bash_variable;
    my $FILEHANDLE;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        bash_variable => {
            defined     => 1,
            required    => 1,
            strict_type => 1,
            store       => \$bash_variable,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path,
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append,
        },
        stdoutfile_path => {
            strict_type => 1,
            store       => \$stdoutfile_path,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores command
    my @commands = q{unset};

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
            FILEHANDLE   => $FILEHANDLE,
            commands_ref => \@commands,
            separator    => $SPACE,

        }
    );
    return @commands;
}
1;
