package MIP::Program::Vt;

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
use autodie qw{ :all };
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
    our @EXPORT_OK = qw{ vt_decompose vt_normalize vt_uniq };
}

Readonly my $BASE_COMMAND => q{vt};

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
            store       => \$infile_path,
            strict_type => 1,
        },
        outfile_path        => { store => \$outfile_path, strict_type => 1, },
        smart_decomposition => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$smart_decomposition,
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
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = (
        get_executable_base_command( { base_command => $BASE_COMMAND, } ),
        qw{ decompose }
    );

    if ($smart_decomposition) {

        push @commands, q{-s};
    }

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
            commands_ref => \@commands,
            filehandle   => $filehandle,
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
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        no_fail_inconsistent_reference => {
            allow       => [ 0, 1 ],
            default     => 0,
            strict_type => 1,
            store       => \$no_fail_inconsistent_reference,
        },
        outfile_path       => { store => \$outfile_path, strict_type => 1, },
        referencefile_path => {
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
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = (
        get_executable_base_command( { base_command => $BASE_COMMAND, } ),
        qw{ normalize }
    );

    if ($no_fail_inconsistent_reference) {

        push @commands, q{-n};
    }

    if ($referencefile_path) {

        push @commands, q{-r} . $SPACE . $referencefile_path;
    }

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
            commands_ref => \@commands,
            filehandle   => $filehandle,
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
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outfile_path    => { store => \$outfile_path, strict_type => 1, },
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        stderrfile_path_append => {
            store       => \$stderrfile_path_append,
            strict_type => 1,
        },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands =
      ( get_executable_base_command( { base_command => $BASE_COMMAND, } ), qw{ uniq } );

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
            commands_ref => \@commands,
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

1;
