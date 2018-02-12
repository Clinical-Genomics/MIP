package MIP::Recipes::Analysis::Xargs;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ fileparse };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ xargs_command };

}

## Constants
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $PIPE       => q{|};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub xargs_command {

##Function : Creates the command line for xargs. Writes to sbatch FILEHANDLE and opens xargs FILEHANDLE
##Returns  : xargs_file_counter, $xargs_file_path
##Arguments: $core_number          => Number of cores to use
##         : $FILEHANDLE           => Sbatch filehandle to write to
##         : $file_path            => File path
##         : $first_command        => Inital command
##         : $java_jar             => Java jar
##         : $java_use_large_pages => Use java large pages {REF}
##         : $memory_allocation    => Memory allocation for java
##         : $program_info_path    => Program info path
##         : $XARGSFILEHANDLE      => XARGS filehandle to write to
##         : $temp_directory       => Redirect tmp files to java temp {Optional}
##         : $xargs_file_counter   => Xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $core_number;
    my $FILEHANDLE;
    my $file_path;
    my $first_command;
    my $java_jar;
    my $java_use_large_pages;
    my $memory_allocation;
    my $program_info_path;
    my $XARGSFILEHANDLE;
    my $temp_directory;

    ## Default(s)
    my $xargs_file_counter;

    my $tmpl = {
        core_number => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$core_number
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE },
        file_path  => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$file_path
        },
        first_command        => { strict_type => 1, store => \$first_command },
        java_jar             => { strict_type => 1, store => \$java_jar },
        java_use_large_pages => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$java_use_large_pages
        },
        memory_allocation => { strict_type => 1, store => \$memory_allocation },
        program_info_path => { strict_type => 1, store => \$program_info_path },
        temp_directory    => { strict_type => 1, store => \$temp_directory },
        XARGSFILEHANDLE =>
          { required => 1, defined => 1, store => \$XARGSFILEHANDLE },
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
            FILEHANDLE       => $FILEHANDLE,
            infile_paths_ref => [$xargs_file_path],
        }
    );

    # Pipe
    print {$FILEHANDLE} $PIPE . $SPACE;

    my @commands;
    if ( ( defined $first_command ) && ( $first_command eq q{java} ) ) {

        ## Return java core commands
        @commands = java_core(
            {
                java_jar             => $java_jar,
                java_use_large_pages => $java_use_large_pages,
                memory_allocation    => $memory_allocation,
                temp_directory       => $temp_directory,
            }
        );
    }
    elsif ($first_command) {

        @commands = ($first_command);
    }

    xargs(
        {
            FILEHANDLE         => $FILEHANDLE,
            max_procs          => $core_number,
            shell_commands_ref => \@commands,
            replace_str        => 1,
            verbose            => 1,
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
