package MIP::Program::Variantcalling::Vt;

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
use Readonly;

## MIPs lib/
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ vt_decompose vt_normalize vt_uniq };
}

## Constants
Readonly my $SPACE => q{ };

sub vt_decompose {

## Function : Perl wrapper for writing Vt decompose recipe to $filehandle or return commands array. Based on Vt v0.5.
## Returns  : @commands
## Arguments: $filehandle             => Filehandle to write to
##          : $infile_path            => Infile path to read from
##          : $outfile_path           => Outfile path to write to
##          : $smart_decomposition    => Smart decomposition
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $smart_decomposition;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path,
        },
        outfile_path        => { strict_type => 1, store => \$outfile_path, },
        smart_decomposition => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$smart_decomposition,
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path,
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append,
        },
        stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ vt decompose };

    if ($smart_decomposition) {

        push @commands, q{-s};
    }

    # Specify output filename
    if ($outfile_path) {

        push @commands, q{-o} . $SPACE . $outfile_path;
    }

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
            filehandle   => $filehandle,
            commands_ref => \@commands,
            separator    => $SPACE,

        }
    );
    return @commands;
}

sub vt_normalize {

## Function : Perl wrapper for writing Vt normalize recipe to $filehandle or return commands array. Based on Vt v0.5.
## Returns  : @commands
## Arguments: $filehandle                     => Filehandle to write to
##          : $infile_path                    => Infile path to read from
##          : $no_fail_inconsistent_reference => Do not fail when REF is inconsistent with reference sequence for non SNPs
##          : $outfile_path                   => Outfile path to write to
##          : $referencefile_path             => Reference sequence fasta file
##          : $stderrfile_path                => Stderrfile path
##          : $stderrfile_path_append         => Append stderr info to file path
##          : $stdoutfile_path                => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $outfile_path;
    my $referencefile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $no_fail_inconsistent_reference;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path,
        },
        no_fail_inconsistent_reference => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$no_fail_inconsistent_reference,
        },
        outfile_path       => { strict_type => 1, store => \$outfile_path },
        referencefile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$referencefile_path,
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path,
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append,
        },
        stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ vt normalize };

    ## Options
    if ($no_fail_inconsistent_reference) {

        push @commands, q{-n};
    }

    if ($referencefile_path) {

        push @commands, q{-r} . $SPACE . $referencefile_path;
    }

    #Specify output filename
    if ($outfile_path) {

        push @commands, q{-o} . $SPACE . $outfile_path;
    }

    ## Infile
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
            filehandle   => $filehandle,
            commands_ref => \@commands,
            separator    => $SPACE,

        }
    );
    return @commands;
}

sub vt_uniq {

## Function : Perl wrapper for writing Vt normalize recipe to $filehandle or return commands array. Based on Vt v0.5.
## Returns  : @commands
## Arguments: $filehandle             => Filehandle to write to
##          : $infile_path            => Infile path to read from
##          : $outfile_path           => Outfile path to write to
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path,
        },
        outfile_path    => { strict_type => 1, store => \$outfile_path, },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path,
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append,
        },
        stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ vt uniq };

    #Specify output filename
    if ($outfile_path) {

        push @commands, q{-o} . $SPACE . $outfile_path;
    }

    ## Infile
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
            filehandle   => $filehandle,
            commands_ref => \@commands,
            separator    => $SPACE,

        }
    );
    return @commands;
}

1;
