package MIP::Gnu::Coreutils;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;    #Allow unicode characters in this script
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use autodie;
use Params::Check qw{ check allow last_error };
$Params::Check::PRESERVE_CASE = 1;    #Do not convert to lower case
use Readonly;

use FindBin qw{ $Bin };                 #Find directory of script
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir };

## MIPs lib/
use lib catdir( dirname($Bin), 'lib' );
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    use base qw{ Exporter };
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.04;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ gnu_cp gnu_rm gnu_mv gnu_mkdir gnu_cat gnu_echo gnu_split gnu_sort gnu_printf gnu_sleep gnu_link };
}

## Constants
Readonly my $SPACE     => q{ };
Readonly my $COMMA     => q{,};
Readonly my $EMPTY_STR => q{};
Readonly my $DOUBLE_QUOTE => q{"};

sub gnu_cp {

##gnu_cp

##Function : Perl wrapper for writing cp recipe to already open $FILEHANDLE or return commands array. Based on cp 8.4
##Returns  : "@commands"
##Arguments: $preserve_attributes_ref, $infile_path, $outfile_path, $stderrfile_path, $FILEHANDLE, $preserve, $recursive, $force, $verbose
##         : $preserve_attributes_ref => Preserve the specified attributes (default:mode,ownership,timestamps), if possible additional attributes: context, links, xattr, all
##         : $infile_path            => Infile path
##         : $outfile_path           => Outfile path
##         : $stderrfile_path        => Stderrfile path
##         : $stderrfile_path_append => Append stderrinfo to file
##         : $FILEHANDLE             => Filehandle to write to
##         : $preserve               => Same as --preserve=mode,ownership,timestamps
##         : $recursive              => Copy directories recursively
##         : $force                  => If an existing destination file cannot be opened, remove it and try again
##         : $verbose                => Verbosity

    my ($arg_href) = @_;

    ## Default(s)
    my $preserve;
    my $recursive;
    my $force;
    my $verbose;

    ## Flatten argument(s)
    my $preserve_attributes_ref;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $FILEHANDLE;

    my $tmpl = {
        preserve_attributes_ref => {
            default     => [],
            strict_type => 1,
            store       => \$preserve_attributes_ref
        },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        outfile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfile_path
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append
        },
        FILEHANDLE => {
            store => \$FILEHANDLE
        },
        recursive => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$recursive
        },
        force => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$force
        },
        preserve => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$preserve
        },
        verbose => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$verbose
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## gnu_cp
    # Stores commands depending on input parameters
    my @commands = q{cp};

    # Preserve the specified attributes
    if ( @{$preserve_attributes_ref} ) {
        push @commands,
          q{--preserve=} . join $COMMA, @{$preserve_attributes_ref};
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

    #Explain what is being done
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
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;
}

sub gnu_mv {

## gnu_mv

##Function : Perl wrapper for writing mv recipe to already open $FILEHANDLE or return commands array. Based on mv 8.4
##Returns  : "@commands"
##Arguments: $infile_path, $outfile_path, $stderrfile_path, $stderrfile_path_append, $FILEHANDLE, $force, $verbose
##         : $infile_path            => Infile path
##         : $outfile_path           => Outfile path
##         : $stderrfile_path        => Stderrfile path
##         : $stderrfile_path_append => Append to stderrinfo to file
##         : $FILEHANDLE             => Filehandle to write to
##         : $force                  => If an existing destination file cannot be opened, remove it and try again
##         : $verbose                => Verbosity

    my ($arg_href) = @_;

    ## Default(s)
    my $force;
    my $verbose;

    ## Flatten argument(s)
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $FILEHANDLE;

    my $tmpl = {
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        outfile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfile_path
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append
        },
        FILEHANDLE => {
            store => \$FILEHANDLE
        },
        force => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$force
        },
        verbose => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$verbose
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## gnu_mv
    # Stores commands depending on input parameters
    my @commands = q{mv};

    if ($force) {
        push @commands, q{--force};
    }

    #Explain what is being done
    if ($verbose) {
        push @commands, q{--verbose};
    }

    push @commands, $infile_path;
    push @commands, $outfile_path;

    #Redirect stderr output to program specific stderr file
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

sub gnu_rm {

##gnu_rm

##Function : Perl wrapper for writing rm recipe to already open $FILEHANDLE or return commands array. Based on rm 8.4
##Returns  : "@commands"
##Arguments: $infile_path, $stderrfile_path, $FILEHANDLE, $recursive, $force, $verbose
##         : $infile_path             => Infile path
##         : $stderrfile_path         => Stderrfile path
##         : $stderrfile_path_append  => Append to stderrinfo to file
##         : $FILEHANDLE              => Filehandle to write to
##         : $recursive               => Copy directories recursively
##         : $force                   => If an existing destination file cannot be opened, remove it and try again
##         : $verbose                 => Verbosity

    my ($arg_href) = @_;

    ## Default(s)
    my $recursive;
    my $force;
    my $verbose;

    ## Flatten argument(s)
    my $infile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $FILEHANDLE;

    my $tmpl = {
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append
        },
        FILEHANDLE => {
            store => \$FILEHANDLE
        },
        recursive => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$recursive
        },
        force => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$force
        },
        verbose => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$verbose
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## gnu_rm
    # Stores commands depending on input parameters
    my @commands = q{rm};

    if ($recursive) {
        push @commands, q{--recursive};
    }

    if ($force) {
        push @commands, q{--force};
    }

    # Explain what is being done
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
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;
}

sub gnu_mkdir {

## gnu_mkdir

##Function : Perl wrapper for writing mkdir recipe to already open $FILEHANDLE or return commands array. Based on mkdir 8.4
##Returns  : "@commands"
##Arguments: $indirectory_path, $stderrfile_path, $stderrfile_path_append, $FILEHANDLE, $parents, $verbose
##         : $indirectory_path        => Infile path
##         : $stderrfile_path         => Stderrfile path
##         : $stderrfile_path_append  => Append to stderrinfo to file
##         : $FILEHANDLE              => Filehandle to write to
##         : $parents                 => No error if existing, make parent directories as needed
##         : $verbose                 => Verbosity

    my ($arg_href) = @_;

    ## Default(s)
    my $verbose;
    my $parents;

    ## Flatten argument(s)
    my $indirectory_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $FILEHANDLE;

    my $tmpl = {
        indirectory_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$indirectory_path
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append
        },
        FILEHANDLE => {
            store => \$FILEHANDLE
        },
        parents => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$parents
        },
        verbose => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$verbose
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## gnu_mkdir
    # Stores commands depending on input parametersr
    my @commands = q{mkdir};

    # Make parent directories as needed
    if ($parents) {
        push @commands, q{--parents};
    }

    # Explain what is being done
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
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;
}

sub gnu_cat {

##gnu_cat

##Function : Perl wrapper for writing cat recipe to already open $FILEHANDLE or return commands array. Based on cat 8.4
##Returns  : "@commands"
##Arguments: $infile_paths_ref, $outfile_path, $stderrfile_path, $stderrfile_path_append, $FILEHANDLE
##         : $infile_paths_ref        => Infile paths {REF}
##         : $outfile_path            => Outfile path
##         : $stderrfile_path         => Stderrfile path
##         : $stderrfile_path_append  => Append to stderrinfo to file
##         : $FILEHANDLE              => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_paths_ref;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $FILEHANDLE;

    my $tmpl = {
        infile_paths_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$infile_paths_ref
        },
        outfile_path => {
            strict_type => 1,
            store       => \$outfile_path
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append
        },
        FILEHANDLE => {
            store => \$FILEHANDLE
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## gnu_cat
    # Stores commands depending on input parameters
    my @commands = q{cat};

    ## Infiles
    push @commands, join $SPACE, @{$infile_paths_ref};

    ## Outfile
    if ($outfile_path) {
        push @commands, q{>} . $SPACE . $outfile_path;
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

sub gnu_echo {

##gnu_echo

##Function : Perl wrapper for writing echo recipe to already open $FILEHANDLE or return commands array. Based on echo 8.4
##Returns  : "@commands"
##Arguments: $strings_ref, $outfile_path, $stderrfile_path,$stderrfile_path_append, $FILEHANDLE, $enable_interpretation, $no_trailing_newline
##         : $strings_ref             => Strings to echo {REF}
##         : $outfile_path            => Outfile path
##         : $stderrfile_path         => Stderrfile path
##         : $stderrfile_path_append  => Append to stderrinfo to file
##         : $FILEHANDLE              => Filehandle to write to
##         : $enable_interpretation   => Enable interpretation of backslash escapes
##         : $no_trailing_newline     => Do not output the trailing newline

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $strings_ref;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $FILEHANDLE;
    my $enable_interpretation;
    my $no_trailing_newline;

    my $tmpl = {
        strings_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$strings_ref
        },
        outfile_path => {
            strict_type => 1,
            store       => \$outfile_path
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append
        },
        FILEHANDLE => {
            store => \$FILEHANDLE
        },
        enable_interpretation => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$enable_interpretation
        },
        no_trailing_newline => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$no_trailing_newline
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Echo
    # Stores commands depending on input parameters
    my @commands = q{echo};

    ##Options
    if ($enable_interpretation) {
        push @commands, q{-e};
    }

    if ($no_trailing_newline) {
        push @commands, q{-n};
    }

    ## Strings
    push @commands, 
      $DOUBLE_QUOTE . join( $EMPTY_STR, @{$strings_ref} ) . $DOUBLE_QUOTE;

    ## Outfile
    if ($outfile_path) {
        push @commands, q{>} . $SPACE . $outfile_path;
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

sub gnu_split {

## gnu_split

##Function : Perl wrapper for writing split recipe to $FILEHANDLE or return commands array. Based on split 8.4.
##Returns  : "@commands"
##Arguments: $infile_path, $FILEHANDLE, $stderrfile_path, $stderrfile_path_append, $prefix, $lines, $suffix_length, $numeric_suffixes, $quiet, $verbose
##         : $infile_path            => Infile path
##         : $FILEHANDLE             => Filehandle to write to
##         : $stderrfile_path        => Stderrfile path
##         : $stderrfile_path_append => Append to stderrinfo to file
##         : $prefix                 => Prefix of output files
##         : $lines                  => Put number lines per output file
##         : $suffix_length          => Use suffixes of length N
##         : $numeric_suffixes       => Use numeric suffixes instead of alphabetic
##         : $quiet                  => Suppress all warnings
##         : $verbose                => Verbosity

    my ($arg_href) = @_;

    ## Default(s)
    my $quiet;
    my $verbose;

    ## Flatten argument(s)
    my $infile_path;
    my $FILEHANDLE;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $prefix;
    my $lines;
    my $suffix_length;
    my $numeric_suffixes;

    my $tmpl = {
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path
        },
        FILEHANDLE => {
            store => \$FILEHANDLE
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append
        },
        prefix => {
            strict_type => 1,
            store       => \$prefix
        },
        lines => {
            allow       => qr/^\d+$/,
            strict_type => 1,
            store       => \$lines
        },
        suffix_length => {
            allow       => qr/^\d+$/,
            strict_type => 1,
            store       => \$suffix_length
        },
        numeric_suffixes => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$numeric_suffixes
        },
        quiet => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$quiet
        },
        verbose => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$verbose
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## gnu_split
    # Stores commands depending on input parameters
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
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;
}

sub gnu_sort {

##gnu_sort

##Function : Perl wrapper for writing sort recipe to already open $FILEHANDLE or return commands array. Based on sort 8.4.
##Returns  : "@commands"
##Arguments: $keys_ref, $infile_path, $outfile_path, $stderrfile_path,$stderrfile_path_append, $stdoutfile_path, $FILEHANDLE
##         : $keys_ref               => Start a key at POS1 (origin 1), end it at POS2
##         : $infile_path            => Infile path
##         : $outfile_path           => Outfile path
##         : $stderrfile_path        => Stderrfile path
##         : $stderrfile_path_append => Append to stderrinfo to file
##         : $stdoutfile_path        => Stdoutfile path
##         : $FILEHANDLE             => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $keys_ref;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $FILEHANDLE;

    my $tmpl = {
        keys_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$keys_ref
        },
        infile_path => {
            strict_type => 1,
            store       => \$infile_path
        },
        outfile_path => {
            strict_type => 1,
            store       => \$outfile_path
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append
        },
        stdoutfile_path => {
            strict_type => 1,
            store       => \$stdoutfile_path
        },
        FILEHANDLE => {
            store => \$FILEHANDLE
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## gnu_sort
    # Stores commands depending on input parameters
    my @commands = q{sort};

    ## Options
    if ( @{$keys_ref} ) {
        push @commands, q{--key} . $SPACE 
          . join $SPACE . q{--key} . $SPACE, @{$keys_ref};
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

sub gnu_printf {

##gnu_printf

##Function : Perl wrapper for writing printf recipe to already open $FILEHANDLE or return commands array. Based on printf 8.4.
##Returns  : "@commands"
##Arguments: $format_string, $stderrfile_path, $stderrfile_path_append, $stdoutfile_path, $FILEHANDLE
##         : $format_string          => Format string to print
##         : $stderrfile_path        => Stderrfile path
##         : $stderrfile_path_append => Append to stderrinfo to file
##         : $stdoutfile_path        => Stdoutfile path
##         : $FILEHANDLE             => Filehandle to write to
##         : $append_stderr_info     => Append stderr info to file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $format_string;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $FILEHANDLE;

    my $tmpl = {
        format_string => {
            strict_type => 1,
            store       => \$format_string
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append
        },
        stdoutfile_path => {
            strict_type => 1,
            store       => \$stdoutfile_path
        },
        FILEHANDLE => {
            store => \$FILEHANDLE
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## gnu_printf
    # Stores commands depending on input parametersf
    my @commands = q{printf};

    ## Options
    if ($format_string) {
        push @commands, $format_string;
    }

    #Redirect stdout to program specific stdout file
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

sub gnu_sleep {

##gnu_sleep

##Function : Perl wrapper for writing sleep recipe to already open $FILEHANDLE or return commands array. Based on sleep 8.4
## Returns : "@commands"
##Arguments: $seconds_to_sleep, $stderrfile_path, $stderrfile_path_append, $stdoutfile_path, $FILEHANDLE
##         : $seconds_to_sleep       => Seconds to sleep
##         : $outfile_path           => Outfile path
##         : $stderrfile_path        => Stderrfile path
##         : $stderrfile_path_append => Append to stderrinfo to file
##         : $stdoutfile_path        => Stdoutfile path
##         : $FILEHANDLE             => Filehandle to write to
##         : $append_stderr_info     => Append stderr info to file

    my ($arg_href) = @_;

    ## Default(s)
    my $seconds_to_sleep;

    ## Flatten argument(s)
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $FILEHANDLE;

    my $tmpl = {
        seconds_to_sleep => {
            default     => 0,
            allow       => qr/^\d+$/,
            strict_type => 1,
            store       => \$seconds_to_sleep
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append
        },
        stdoutfile_path => {
            strict_type => 1,
            store       => \$stdoutfile_path
        },
        FILEHANDLE => {
            store => \$FILEHANDLE
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## gnu_sleep
    # Stores commands depending on input parametersf
    my @commands = q{sleep};

    ## Options
    if ($seconds_to_sleep) {

        push @commands, $seconds_to_sleep;
    }

    #Redirect stdout to program specific stdout file
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

sub gnu_link {

##gnu_ln

##Function : Perl wrapper for writing link recipe to already open $FILEHANDLE or return commands array. Based on link 8.4
## Returns : "@commands"
##Arguments: $symbolic, $force, $link_path, $target_path, $stderrfile_path, $stderrfile_path_append, $stdoutfile_path, $FILEHANDLE
##         : $symbolic               => Create a symbolic link
##         : $force                  => Remove existing destination files
##         : $link_path              => Path to link
##         : $target_path            => Path to target
##         : $stderrfile_path        => Stderrfile path
##         : $stderrfile_path_append => Append to stderrinfo to file
##         : $stdoutfile_path        => Stdoutfile path
##         : $FILEHANDLE             => Filehandle to write to

    my ($arg_href) = @_;

    ## Default(s)

    ## Flatten argument(s)
    my $symbolic;
    my $force;
    my $link_path;
    my $target_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $FILEHANDLE;

    my $tmpl = {
        symbolic => {
            default => 0,
            allow => [0, 1],
            strict_type => 1,
            store => \$symbolic
        },
        force => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$force
        },
        link_path => {
            required => 1,
            strict_type => 1,
            store => \$link_path
        },
        target_path => {
            required => 1,
            strict_type => 1,
            store => \$target_path,
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append
        },
        stdoutfile_path => {
            strict_type => 1,
            store       => \$stdoutfile_path
        },
        FILEHANDLE => {
            store => \$FILEHANDLE
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## gnu_link
    # Stores commands depending on input parametersf
    my @commands = q{link};

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
