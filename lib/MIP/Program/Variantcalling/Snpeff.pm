package MIP::Program::Variantcalling::Snpeff;

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
use MIP::Gnu::Coreutils qw{ gnu_rm };
use MIP::Language::Java qw{ java_core };
use MIP::Script::Utils qw{ create_temp_dir };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ snpeff_ann snpeff_download };
}

## Constants
Readonly my $DASH    => q{-};
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };

sub snpeff_ann {

## Function : Perl wrapper for writing snpeff ann recipe to already open $FILEHANDLE or return commands array. Based on SnpEff 4.2 (build 2015-12-05).
## Returns  : @commands
## Arguments: $config_file_path       => Config file path
##          : $FILEHANDLE             => Filehandle to write to
##          : $genome_build_version   => Genome build version
##          : $infile_path            => Infile path
##          : $outfile_path           => Outfile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $verbosity              => Increase output verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $config_file_path;
    my $FILEHANDLE;
    my $genome_build_version;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $verbosity;

    my $tmpl = {
        config_file_path => { strict_type => 1, store => \$config_file_path },
        FILEHANDLE       => {
            store => \$FILEHANDLE,
        },
        genome_build_version => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$genome_build_version,
        },
        infile_path     => { strict_type => 1, store => \$infile_path, },
        outfile_path    => { strict_type => 1, store => \$outfile_path, },
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
        verbosity => {
            allow       => qr/^\w+$/,
            strict_type => 1,
            store       => \$verbosity,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = q{ann};

    ## Options
    if ($verbosity) {

        push @commands, $DASH . $verbosity;
    }

    if ($genome_build_version) {

        push @commands, $genome_build_version;
    }

    if ($config_file_path) {

        push @commands, q{-config} . $SPACE . $config_file_path;
    }

    ## Infile
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
            FILEHANDLE   => $FILEHANDLE,
            commands_ref => \@commands,
            separator    => $SPACE,

        }
    );
    return @commands;
}

sub snpeff_download {

## Function : Write instructions to download snpeff database
## Returns  : @commands
## Arguments: $config_file_path        => Path to snpeff config file
##          : $FILEHANDLE              => FILEHANDLE to write to
##          : $genome_version_database => Database to be downloaded
##          : $jar_path                => Path to snpeff jar
##          : $memory_allocation       => Java memory allocation
##          : $stderrfile_path         => Stderrfile path
##          : $stderrfile_path_append  => Append to stderrinfo to file
##          : $stdoutfile_path         => Stdoutfile path
##          : $temp_directory          => Temporary directory
##          : $verbose                 => Verbose output

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $config_file_path;
    my $FILEHANDLE;
    my $genome_version_database;
    my $jar_path;
    my $memory_allocation;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $temp_directory;
    my $verbose;

    my $tmpl = {
        config_file_path => {
            defined     => 1,
            strict_type => 1,
            store       => \$config_file_path,
        },
        FILEHANDLE => {
            required => 1,
            store    => \$FILEHANDLE,
        },
        genome_version_database => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$genome_version_database,
        },
        jar_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$jar_path,
        },
        memory_allocation => {
            default     => q{Xmx2g},
            defined     => 1,
            strict_type => 1,
            store       => \$memory_allocation,
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
        temp_directory => {
            default => 0,
            allow   => [ 0, 1 ],
            store   => \$temp_directory,
        },
        verbose => {
            default     => 1,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$verbose,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Create optional temporary directory
    if ($temp_directory) {
        $temp_directory = create_temp_dir( { FILEHANDLE => $FILEHANDLE } );
        say {$FILEHANDLE} $NEWLINE;
    }

    ## Build base command
    my @base = java_core(
        {
            memory_allocation => $memory_allocation,
            java_jar          => $jar_path,
            temp_directory    => $temp_directory,
        }
    );
    my @commands = ( @base, qw{download} );

    ## Add database to be downloaded
    push @commands, $genome_version_database;

    ## Add verbose flag
    if ($verbose) {
        push @commands, q{-v};
    }

    ## Add otional path to config file
    if ($config_file_path) {
        push @commands, q{-c} . $SPACE . $config_file_path;
    }

    ## Optionally capture output
    push @commands,
      unix_standard_streams(
        {
            stdoutfile_path        => $stdoutfile_path,
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
        }
      );

    ## Write rest of java commadn to $FILEHANDLE
    unix_write_to_file(
        {
            commands_ref => \@commands,
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    if ($temp_directory) {
        say {$FILEHANDLE} $NEWLINE;
        gnu_rm(
            {
                infile_path => $temp_directory,
                recursive   => 1,
                FILEHANDLE  => $FILEHANDLE,
            }
        );
    }

    return @commands;
}

1;
