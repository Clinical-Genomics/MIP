package MIP::Program::Cyrius;

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
use MIP::Constants qw{ $SPACE };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

Readonly my $GRCH37 => 37;
Readonly my $GRCH38 => 38;
Readonly my $HG19   => 19;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ star_caller };
}

sub star_caller {

## Function : Perl wrapper for star_caller.py. Based on Cyrius v1.0
## Returns  : @commands
## Arguments: $filehandle             => Filehandle to write to
##          : $include_new_star       => Include newly added uncurated star alleles
##          : $known_function         => Only call alleles with known function
##          : $manifest_file_path     => List of sample infiles
##          : $outdir_path            => Out directory path
##          : $outfile_prefix         => File prefix
##          : $reference_version      => Genome reference version
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdinfile_path         => Stdinfile path
##          : $stdoutfile_path        => Stdoutfile path
##          : $thread_number          => Threads

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $genome_version;
    my $include_new_star;
    my $known_function;
    my $manifest_file_path;
    my $outdir_path;
    my $outfile_prefix;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdinfile_path;
    my $stdoutfile_path;
    my $thread_number;

    ## Default(s)

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        genome_version => {
            allow       => [ $GRCH37, $GRCH38, $HG19 ],
            required    => 1,
            store       => \$genome_version,
            strict_type => 1,
        },
        include_new_star => {
            allow       => [ undef, 0, 1 ],
            store       => \$include_new_star,
            strict_type => 1,
        },
        known_function => {
            allow       => [ undef, 0, 1 ],
            store       => \$known_function,
            strict_type => 1,
        },
        manifest_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$manifest_file_path,
            strict_type => 1,
        },
        outdir_path => {
            defined     => 1,
            required    => 1,
            store       => \$outdir_path,
            strict_type => 1,
        },
        outfile_prefix => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_prefix,
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
        thread_number => {
            allow       => qr/\A \d+ \z/xms,
            store       => \$thread_number,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ star_caller.py };

    push @commands, q{--genome} . $SPACE . $genome_version;

    if ($include_new_star) {

        push @commands, q{--includeNewStar};
    }

    if ($known_function) {

        push @commands, q{--knownFunction};
    }

    push @commands, q{--manifest} . $SPACE . $manifest_file_path;

    push @commands, q{--outDir} . $SPACE . $outdir_path;

    push @commands, q{--prefix} . $SPACE . $outfile_prefix;

    if ($thread_number) {

        push @commands, q{--threads} . $SPACE . $thread_number;
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
