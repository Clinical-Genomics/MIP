package MIP::Gnu::Coreutils;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use FindBin qw{ $Bin };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;

## CPANM
use autodie;
use Readonly;

## MIPs lib/
use MIP::Constants
  qw{ $AMPERSAND $COMMA $DOUBLE_QUOTE $EMPTY_STR $EQUALS $NEWLINE $SINGLE_QUOTE $SPACE };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    use base qw{ Exporter };
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.12;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ gnu_cat
      gnu_chmod
      gnu_cp
      gnu_cut
      gnu_echo
      gnu_head
      gnu_ln
      gnu_md5sum
      gnu_mkdir
      gnu_mv
      gnu_printf
      gnu_rm
      gnu_rm_and_echo
      gnu_sleep
      gnu_sort
      gnu_split
      gnu_tail
      gnu_touch
      gnu_uniq
    };
}

sub gnu_cat {

## Function : Perl wrapper for writing cat command to already open $FILEHANDLE or return commands array. Based on cat 8.4
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $infile_paths_ref       => Infile paths {REF}
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append to stderrinfo to file
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_paths_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        infile_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$infile_paths_ref,
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
    my @commands = qw{ cat };

    ## Infiles
    push @commands, join $SPACE, @{$infile_paths_ref};

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

sub gnu_chmod {

## Function : Perl wrapper for writing chmod command to already open $FILEHANDLE or return commands array. Based on chmod 8.4
## Returns  : @commands
## Arguments: $FILEHANDLE             => FILEHANDLE to write to
##          : $file_path              => Path to file
##          : $permission             => Permisions for the file
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append to stderrinfo to file
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $file_path;
    my $permission;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        permission => {
            defined     => 1,
            required    => 1,
            store       => \$permission,
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

    my @commands = qw{ chmod };

    ## Add the permission
    push @commands, $permission;

    ## Add the file path
    push @commands, $file_path;

    ## Redirect stdout to program specific stdout file
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

sub gnu_cp {

## Function : Perl wrapper for writing cp command to already open $FILEHANDLE or return commands array. Based on cp 8.4
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $force                  => If an existing destination file cannot be opened, remove it and try again
##          : $infile_path            => Infile path
##          : $outfile_path           => Outfile path
##          : $preserve               => Same as --preserve=mode,ownership,timestamps
##          : $recursive              => Copy directories recursively
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderrinfo to file
##          : $verbose                => Verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_path;
    my $preserve_attributes_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;

    ## Default(s)
    my $force;
    my $preserve;
    my $recursive;
    my $verbose;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        force => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$force,
            strict_type => 1,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        preserve => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$preserve,
            strict_type => 1,
        },
        preserve_attributes_ref => {
            default     => [],
            store       => \$preserve_attributes_ref,
            strict_type => 1,
        },
        recursive => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$recursive,
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
        verbose => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$verbose,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ cp };

    ## Preserve the specified attributes
    if ( @{$preserve_attributes_ref} ) {
        push @commands, q{--preserve=} . join $COMMA, @{$preserve_attributes_ref};
    }

    elsif ($preserve) {
        push @commands, q{-p};
    }

    if ($recursive) {
        push @commands, q{--recursive};
    }

    if ($force) {
        push @commands, q{--force};
    }

    ## Explain what is being done
    if ($verbose) {
        push @commands, q{--verbose};
    }

    push @commands, $infile_path;
    push @commands, $outfile_path;

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

sub gnu_cut {

## Function : Perl wrapper for writing cut command to already open $FILEHANDLE or return commands array. Based on cut 8.4
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $infile_path            => Infile paths {REF}
##          : $list                   => List of specified fields
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append to stderrinfo to file
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $list;
    my $infile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        list => {
            store       => \$list,
            strict_type => 1,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
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
    my @commands = qw{ cut };

    if ($list) {

        push @commands, q{-f} . $SPACE . $list;
    }

    ## Infiles
    push @commands, $infile_path;

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

sub gnu_echo {

## Function : Perl wrapper for writing echo command to already open $FILEHANDLE or return commands array. Based on echo 8.4
## Returns  : @commands
## Arguments: $enable_interpretation  => Enable interpretation of backslash escapes
##          : $FILEHANDLE             => Filehandle to write to
##          : $no_trailing_newline    => Do not output the trailing newline
##          : $outfile_path           => Outfile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append to stderrinfo to file
##          : $stdoutfile_path_append => Append to stdout to file
##          : $strings_ref            => Strings to echo {REF}
##          : $string_wrapper         => Wrap string with this character

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $enable_interpretation;
    my $FILEHANDLE;
    my $no_trailing_newline;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path_append;
    my $strings_ref;
    my $string_wrapper;

    my $tmpl = {
        enable_interpretation => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$enable_interpretation,
            strict_type => 1,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        no_trailing_newline => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$no_trailing_newline,
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
        stdoutfile_path_append => {
            store       => \$stdoutfile_path_append,
            strict_type => 1,
        },
        strings_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$strings_ref,
            strict_type => 1,
        },
        string_wrapper => {
            allow       => [ $DOUBLE_QUOTE, $SINGLE_QUOTE ],
            default     => $DOUBLE_QUOTE,
            defined     => 1,
            store       => \$string_wrapper,
            strict_type => 1,
        }

    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = q{echo};

    ## Options
    if ($enable_interpretation) {
        push @commands, q{-e};
    }

    if ($no_trailing_newline) {
        push @commands, q{-n};
    }

    ## Strings
    push @commands,
      $string_wrapper . join( $EMPTY_STR, @{$strings_ref} ) . $string_wrapper;

    ## Outfile
    if ($outfile_path) {
        push @commands, q{>} . $SPACE . $outfile_path;
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
            stdoutfile_path_append => $stdoutfile_path_append,
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

sub gnu_head {

## Function : Perl wrapper for writing head command to already open $FILEHANDLE or return commands array
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $infile_path            => Infile path
##          : $lines                  => Lines to print
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $lines;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        infile_path => {
            store       => \$infile_path,
            strict_type => 1,
        },
        lines => {
            allow       => qr/ ^\d+$ /xms,
            store       => \$lines,
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
    my @commands = q{head};

    if ($lines) {
        push @commands, q{-n} . $SPACE . $lines;
    }

    if ($infile_path) {

        push @commands, $infile_path;
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

sub gnu_ln {

## Function : Perl wrapper for writing ln command to already open $FILEHANDLE or return commands array. Based on ln 8.4
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $force                  => Remove existing destination files
##          : $link_path              => Path to link
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append to stderrinfo to file
##          : $stdoutfile_path        => Stdoutfile path
##          : $symbolic               => Create a symbolic link
##          : $target_path            => Path to target

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $force;
    my $link_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $symbolic;
    my $target_path;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        force => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$force,
            strict_type => 1,
        },
        link_path => {
            required    => 1,
            store       => \$link_path,
            strict_type => 1,
        },
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append,
        },
        stdoutfile_path => {
            store       => \$stdoutfile_path,
            strict_type => 1,
        },
        symbolic => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$symbolic,
            strict_type => 1,
        },
        target_path => {
            required    => 1,
            store       => \$target_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parametersf
    my @commands = q{ln};

    ## Options
    if ($symbolic) {
        push @commands, q{--symbolic};
    }

    if ($force) {
        push @commands, q{--force};
    }

    #Add target and link path
    push @commands, $target_path;

    push @commands, $link_path;

    #Redirect stdout to program specific stdout file
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

sub gnu_md5sum {

## Function : Perl wrapper for writing md5sum command to already open $FILEHANDLE or return commands array. Based on md5sum 8.4
## Returns  : @commands
## Arguments: $check                  => Read MD5 sums from the FILEs and check them
##          : $FILEHANDLE             => Filehandle to write to
##          : $infile_path            => Infile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $check;
    my $FILEHANDLE;
    my $infile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)

    my $tmpl = {
        check => {
            allow       => [ undef, 0, 1 ],
            store       => \$check,
            strict_type => 1,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        infile_path => {
            store       => \$infile_path,
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
    my @commands = q{md5sum};

    if ($check) {

        push @commands, q{--check};
    }

    if ($infile_path) {

        push @commands, $infile_path;
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

sub gnu_mkdir {

## Function : Perl wrapper for writing mkdir command to already open $FILEHANDLE or return commands array. Based on mkdir 8.4
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $indirectory_path       => Infile path
##          : $parents                => No error if existing, make parent directories as needed
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append to stderrinfo to file
##          : $verbose                => Verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $indirectory_path;
    my $stderrfile_path;
    my $stderrfile_path_append;

    ## Default(s)
    my $parents;
    my $verbose;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        indirectory_path => {
            defined     => 1,
            required    => 1,
            store       => \$indirectory_path,
            strict_type => 1,
        },
        parents => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$parents,
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
        verbose => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$verbose,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parametersr
    my @commands = q{mkdir};

    ## Make parent directories as needed
    if ($parents) {
        push @commands, q{--parents};
    }

    ## Explain what is being done
    if ($verbose) {
        push @commands, q{--verbose};
    }

    ## Indirectory
    push @commands, $indirectory_path;

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

sub gnu_mv {

## Function : Perl wrapper for writing mv command to already open $FILEHANDLE or return commands array. Based on mv 8.4
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $force                  => If an existing destination file cannot be opened, remove it and try again
##          : $infile_path            => Infile path
##          : $outfile_path           => Outfile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append to stderrinfo to file
##          : $verbose                => Verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;

    ## Default(s)
    my $force;
    my $verbose;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        force => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$force,
            strict_type => 1,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
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
        verbose => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$verbose,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = q{mv};

    if ($force) {
        push @commands, q{--force};
    }

    ## Explain what is being done
    if ($verbose) {
        push @commands, q{--verbose};
    }

    push @commands, $infile_path;

    push @commands, $outfile_path;

    ## Redirect stderr output to program specific stderr file
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

sub gnu_printf {

## Function : Perl wrapper for writing printf command to already open $FILEHANDLE or return commands array. Based on printf 8.4.
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $format_string          => Format string to print
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append to stderrinfo to file
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $format_string;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        format_string => {
            store       => \$format_string,
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

    ## Stores commands depending on input parametersf
    my @commands = q{printf};

    ## Options
    if ($format_string) {
        push @commands, $format_string;
    }

    #Redirect stdout to program specific stdout file
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

sub gnu_rm {

## Function : Perl wrapper for writing rm command to already open $FILEHANDLE or return commands array. Based on rm 8.4
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $force                  => If an existing destination file cannot be opened, remove it and try again
##          : $infile_path            => Infile path
##          : $recursive              => Remove directories recursively
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append to stderrinfo to file
##          : $verbose                => Verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;

    ## Default(s)
    my $force;
    my $recursive;
    my $verbose;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        force => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$force,
            strict_type => 1,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
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
        recursive => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$recursive,
            strict_type => 1,
        },
        verbose => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$verbose,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = q{rm};

    if ($recursive) {
        push @commands, q{--recursive};
    }

    if ($force) {
        push @commands, q{--force};
    }

    ## Explain what is being done
    if ($verbose) {
        push @commands, q{--verbose};
    }

    ## Infile
    push @commands, $infile_path;

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

sub gnu_sleep {

## Function : Perl wrapper for writing sleep command to already open $FILEHANDLE or return commands array. Based on sleep 8.4
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $seconds_to_sleep       => Seconds to sleep
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append to stderrinfo to file
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $FILEHANDLE;

    ## Default(s)
    my $seconds_to_sleep;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        seconds_to_sleep => {
            allow       => qr/ ^\d+$ /xms,
            default     => 0,
            store       => \$seconds_to_sleep,
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

    # Stores commands depending on input parameters
    my @commands = q{sleep};

    ## Options
    if ( defined $seconds_to_sleep ) {

        push @commands, $seconds_to_sleep;
    }

    #Redirect stdout to program specific stdout file
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

sub gnu_sort {

## Function : Perl wrapper for writing sort command to already open $FILEHANDLE or return commands array. Based on sort 8.4.
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $infile_path            => Infile path
##          : $keys_ref               => Start a key at POS1 (origin 1), end it at POS2
##          : $outfile_path           => Outfile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append to stderrinfo to file
##          : $stdoutfile_path        => Stdoutfile path
##          : $temporary_directory    => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $keys_ref;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $temporary_directory;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        infile_path => {
            store       => \$infile_path,
            strict_type => 1,
        },
        keys_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$keys_ref,
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
        temporary_directory => {
            store       => \$temporary_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = q{sort};

    ## Options
    if ( @{$keys_ref} ) {
        push @commands, q{--key} . $SPACE . join $SPACE . q{--key} . $SPACE, @{$keys_ref};
    }

    if ($temporary_directory) {

        push @commands, q{--temporary-directory} . $EQUALS . $temporary_directory;
    }
    ## Infile
    if ($infile_path) {
        push @commands, $infile_path;
    }

    ## Outfile
    if ($outfile_path) {
        push @commands, q{>} . $SPACE . $outfile_path;
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

sub gnu_split {

## Function : Perl wrapper for writing split command to $FILEHANDLE or return commands array. Based on split 8.4.
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $infile_path            => Infile path
##          : $lines                  => Put number lines per output file
##          : $numeric_suffixes       => Use numeric suffixes instead of alphabetic
##          : $prefix                 => Prefix of output files
##          : $quiet                  => Suppress all warnings
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append to stderrinfo to file
##          : $suffix_length          => Use suffixes of length N
##          : $verbose                => Verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $lines;
    my $numeric_suffixes;
    my $prefix;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $suffix_length;

    ## Default(s)
    my $quiet;
    my $verbose;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        lines => {
            allow       => qr/ ^\d+$ /xms,
            store       => \$lines,
            strict_type => 1,
        },
        numeric_suffixes => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$numeric_suffixes,
            strict_type => 1,
        },
        prefix => {
            store       => \$prefix,
            strict_type => 1,
        },
        quiet => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$quiet,
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
        suffix_length => {
            allow       => qr/ ^\d+$ /xms,
            store       => \$suffix_length,
            strict_type => 1,
        },
        verbose => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$verbose,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = q{split};

    ## Options
    if ($lines) {
        push @commands, q{--lines=} . $lines;
    }

    if ($numeric_suffixes) {
        push @commands, q{--numeric-suffixes};
    }

    if ($suffix_length) {
        push @commands, q{--suffix-length=} . $suffix_length;
    }

    if ($quiet) {
        push @commands, q{--quiet};
    }

    if ($verbose) {
        push @commands, q{--verbose};
    }

    ## Infile
    push @commands, $infile_path;

    if ($prefix) {
        push @commands, $prefix;
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

sub gnu_tail {

## Function : Perl wrapper for writing tail command to already open $FILEHANDLE or return commands array. Based on tail 8.4
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $lines                  => Lines to print
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $lines;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        lines => {
            allow       => qr/ ^\d+$ /xms,
            store       => \$lines,
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
    my @commands = q{tail};

    if ($lines) {
        push @commands, q{--lines=} . $lines;
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

sub gnu_touch {

## Function : Perl wrapper for writing touch command to already open $FILEHANDLE or return commands array. Based on touch 8.22
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $file                   => File to touch
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $file;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        file => {
            required    => 1,
            store       => \$file,
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
    my @commands = qw{ touch };

    if ($file) {
        push @commands, $file;
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

sub gnu_uniq {

## Function : Perl wrapper for writing uniq command to already open $FILEHANDLE or return commands array. Based on uniq 8.4
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $infile_path            => Infile paths {REF}
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append to stderrinfo to file
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
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
    my @commands = qw{ uniq };

    ## Infiles
    push @commands, $infile_path;

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

sub gnu_rm_and_echo {

## Function : Perl wrapper for writing rm and echo command to already open $FILEHANDLE;
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : file_href               => Hash with file paths and new content {REF}
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append to stderrinfo to file
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $file_href;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $force;

    my $tmpl = {
        FILEHANDLE => {
            required => 1,
            store    => \$FILEHANDLE,
        },
        file_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_href,
            strict_type => 1,
        },
        force => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$force,
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

  FILE_PATH:
    foreach my $file_path ( keys %{$file_href} ) {

        my @commands = gnu_rm(
            {
                infile_path => $file_path,
                force       => $force,
            }
        );

        push @commands, $AMPERSAND . $AMPERSAND;

        push @commands,
          gnu_echo(
            {
                outfile_path => $file_path,
                strings_ref  => [ $file_href->{$file_path} ],
            }
          );

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
        print {$FILEHANDLE} $NEWLINE;
    }

    return 1;
}

1;
