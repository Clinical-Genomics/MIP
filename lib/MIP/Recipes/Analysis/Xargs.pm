package MIP::Recipes::Analysis::Xargs;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use autodie qw{ :all };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use File::Basename qw{ fileparse };
use File::Spec::Functions qw{ catdir catfile };

## CPANM
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ xargs_command };

}

##Constants
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $PIPE       => q{|};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub xargs_command {

##xargs_command

##Function : Creates the command line for xargs. Writes to sbatch FILEHANDLE and opens xargs FILEHANDLE
##Returns  : "xargs_file_counter, $xargs_file_path"
##Arguments: $FILEHANDLE, $XARGSFILEHANDLE, $file_path, $core_number, $first_command, $program_info_path, $memory_allocation, $java_use_large_pages, $temp_directory, $java_jar, $xargs_file_counter
##         : $FILEHANDLE           => Sbatch filehandle to write to
##         : $XARGSFILEHANDLE      => XARGS filehandle to write to
##         : $file_path            => File path
##         : $core_number          => Number of cores to use
##         : $first_command        => Inital command
##         : $program_info_path    => Program info path
##         : $memory_allocation    => Memory allocation for java
##         : $java_use_large_pages => Use java large pages {REF}
##         : $temp_directory       => Redirect tmp files to java temp {Optional}
##         : $java_jar             => Java jar
##         : $xargs_file_counter   => Xargs file counter

    my ($arg_href) = @_;

    ## Default(s)
    my $xargs_file_counter;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $XARGSFILEHANDLE;
    my $file_path;
    my $core_number;
    my $first_command;
    my $program_info_path;
    my $memory_allocation;
    my $java_use_large_pages;
    my $temp_directory;
    my $java_jar;

    my $tmpl = {
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
        XARGSFILEHANDLE =>
          { required => 1, defined => 1, store => \$XARGSFILEHANDLE },
        file_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$file_path
        },
        core_number => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$core_number
        },
        first_command     => { strict_type => 1, store => \$first_command },
        program_info_path => { strict_type => 1, store => \$program_info_path },
        memory_allocation => { strict_type => 1, store => \$memory_allocation },
        java_use_large_pages => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$java_use_large_pages
        },
        temp_directory     => { strict_type => 1, store => \$temp_directory },
        java_jar           => { strict_type => 1, store => \$java_jar },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$xargs_file_counter
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Gnu::Coreutils qw(gnu_cat);
    use MIP::Gnu::Findutils qw(xargs);
    use MIP::Language::Java qw{java_core};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## File
    my $xargs_file_number = $DOT . $xargs_file_counter;
    my $xargs_file_suffix = $DOT . q{xargs};

    ## Path
    my $xargs_file_path = $file_path . $xargs_file_number . $xargs_file_suffix;

    ## Path
    my $xargs_file_path_prefix;

    ## Check if there is a xargs_file_path_prefix to concatenate
    if ( defined $program_info_path ) {

        $xargs_file_path_prefix = $program_info_path . $xargs_file_number;
    }

    ## Read xargs command file
    gnu_cat(
        {
            infile_paths_ref => [$xargs_file_path],
            FILEHANDLE       => $FILEHANDLE,
        }
    );

    # Pipe
    print {$FILEHANDLE} $PIPE . $SPACE;

    my @commands;
    if ( ( defined $first_command ) && ( $first_command eq q{java} ) ) {

        ## Return java core commands
        @commands = java_core(
            {
                memory_allocation    => $memory_allocation,
                java_use_large_pages => $java_use_large_pages,
                temp_directory       => $temp_directory,
                java_jar             => $java_jar,
            }
        );
    }
    elsif ($first_command) {

        @commands = ($first_command);
    }

    xargs(
        {
            shell_commands_ref => \@commands,
            replace_str        => 1,
            verbose            => 1,
            max_procs          => $core_number,
            FILEHANDLE         => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    # Open xargs file for writing
    open $XARGSFILEHANDLE, q{>}, $xargs_file_path
      or $log->logdie( q{Cannot write to '}
          . $xargs_file_path . q{' :}
          . $OS_ERROR
          . $NEWLINE x 2 );

# Increment to not overwrite xargs file with next call (if used) and $xargs_file_path_prefix
    $xargs_file_counter++;
    return ( $xargs_file_counter, $xargs_file_path_prefix );
}

1;
