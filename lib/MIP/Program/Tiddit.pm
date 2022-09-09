package MIP::Program::Tiddit;

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
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $SPACE };
use MIP::Environment::Executable qw{ get_executable_base_command };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ tiddit_coverage tiddit_sv };

}

## Constants
Readonly my $BASE_COMMAND => q{tiddit};
Readonly my $BIN_SIZE     => 500;

sub tiddit_coverage {

## Function : Perl wrapper for Tiddit coverage. Based on Tiddit 2.7.1.
## Returns  : @commands
## Arguments: $bin_size               => Bin size
##          : $filehandle             => Filehandle to write to
##          : $infile_path            => Infile path
##          : $outfile_path_prefix    => Outfile path prefix
##          : $output_wig             => Generate wig instead of bed
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $outfile_path_prefix;
    my $skip_quality_track;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $bin_size;
    my $output_wig;

    my $tmpl = {
        bin_size => {
            allow       => [ undef, qr/ \A \d+ \z /sxm ],
            default     => $BIN_SIZE,
            store       => \$bin_size,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outfile_path_prefix => { store => \$outfile_path_prefix, strict_type => 1, },
        output_wig          => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$output_wig,
            strict_type => 1,
        },
        skip_quality_track => {
            allow       => [ undef, 0, 1 ],
            store       => \$skip_quality_track,
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

    my @commands =
      ( get_executable_base_command( { base_command => $BASE_COMMAND, } ), qw{ --cov } );

    if ($outfile_path_prefix) {

        push @commands, q{-o} . $SPACE . $outfile_path_prefix;
    }

    if ($bin_size) {

        push @commands, q{-z} . $SPACE . $bin_size;
    }

    if ($output_wig) {

        push @commands, q{-w};
    }

    push @commands, q{--bam} . $SPACE . $infile_path;

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

sub tiddit_sv {

## Function : Perl wrapper for writing tiddit sv recipe to $filehandle or return commands array. Based on tiddit 2.7.1.
## Returns  : @commands
## Arguments: $filehandle                      => Filehandle to write to
##          : $infile_path                     => Infile path
##          : $minimum_number_supporting_pairs => Minimum number of supporting pairs in order to call a variation event
##          : $outfile_path_prefix             => Outfile path. Write documents to FILE
##          : $referencefile_path              => Genome reference file
##          : $stderrfile_path                 => Stderrfile path
##          : $stderrfile_path_append          => Append stderr info to file path
##          : $stdoutfile_path                 => Stdoutfile path
##          : $threads                         => Threads

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $minimum_number_supporting_pairs;
    my $outfile_path_prefix;
    my $referencefile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $threads;

    my $tmpl = {
        filehandle  => { store => \$filehandle, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        minimum_number_supporting_pairs => {
            allow       => qr/ \A \d+ \z /sxm,
            store       => \$minimum_number_supporting_pairs,
            strict_type => 1,
        },
        outfile_path_prefix => { store => \$outfile_path_prefix, strict_type => 1, },
        referencefile_path  => {
            defined     => 1,
            required    => 1,
            store       => \$referencefile_path,
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
        threads => {
            allow       => qr/ \A \d+ \z /xms,
            store       => \$threads,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands =
      ( get_executable_base_command( { base_command => $BASE_COMMAND, } ), qw{ --sv } );

    if ($minimum_number_supporting_pairs) {

        push @commands, q{-p} . $SPACE . $minimum_number_supporting_pairs;
    }

    if ($outfile_path_prefix) {

        push @commands, q{-o} . $SPACE . $outfile_path_prefix;
    }

    push @commands, q{--ref} . $SPACE . $referencefile_path;

    push @commands, q{--bam} . $SPACE . $infile_path;

    if ($threads) {
        push @commands, q{--threads} . $SPACE . $threads;
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
            stdoutfile_path        => $stdoutfile_path,
        }
      );

    push @commands,
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
