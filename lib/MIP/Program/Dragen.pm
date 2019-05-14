package MIP::Program::Dragen;

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

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ dragen_build_hash_table };
}

sub dragen_build_hash_table {

## Function : Perl wrapper for a dragen builing a dragen hash table. Dragen version 3.3.5.
## Returns  : @commands
## Arguments: $build_hash_table           => Build hash table
##          : $FILEHANDLE                 => Filehandle to write to
##          : $ht_alt_liftover_file_path  => Path to lift over file
##          : $ht_decoys_file_path        => Path to decoys file
##          : $outdirectory_path          => Outdirectory path
##          : $reference_genome_file_path => Reference genome file path [.fasta]
##          : $stderrfile_path            => Stderrfile path
##          : $stderrfile_path_append     => Append stderr info to file path
##          : $stdinfile_path             => Stdinfile path
##          : $stdoutfile_path            => Stdoutfile path
##          : $thread_number              => Number of threads

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $build_hash_table;
    my $FILEHANDLE;
    my $ht_alt_liftover_file_path;
    my $ht_decoys_file_path;
    my $outdirectory_path;
    my $reference_genome_file_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdinfile_path;
    my $stdoutfile_path;
    my $thread_number;

    ## Default(s)

    my $tmpl = {
        build_hash_table => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$build_hash_table,
            strict_type => 1,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        ht_alt_liftover_file_path => {
            store       => \$ht_alt_liftover_file_path,
            strict_type => 1,
        },
        ht_decoys_file_path => {
            store       => \$ht_decoys_file_path,
            strict_type => 1,
        },
        outdirectory_path => {
            defined     => 1,
            required    => 1,
            store       => \$outdirectory_path,
            strict_type => 1,
        },
        reference_genome_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$reference_genome_file_path,
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
            allow       => qr{ \A\d+\z }sxm,
            default     => 32,
            store       => \$thread_number,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ dragen };

    return if ( not $build_hash_table );

    push @commands, q{--build-hash-table} . $SPACE . q{true};

    push @commands, q{--ht-num-threads} . $SPACE . $thread_number;

    ## Should match $thread_number
    push @commands, q{--ht-max-table-chunks} . $SPACE . $thread_number;

    push @commands, q{--ht-reference} . $SPACE . $reference_genome_file_path;

    if ($ht_alt_liftover_file_path) {

        push @commands, q{--ht-alt-liftover} . $SPACE . $ht_alt_liftover_file_path;
    }
    if ($ht_decoys_file_path) {

        push @commands, q{--ht-decoys} . $SPACE . $ht_decoys_file_path;
    }
    push @commands, q{--output-directory} . $SPACE . $outdirectory_path;

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
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,

        }
    );
    return @commands;
}

1;
