package MIP::Program::Pdfmerger;

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

    # Set the version for version checking
    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ pdfmerger };
}

Readonly my $BASE_COMMAND => q{pdfmerger};

sub pdfmerger {

## Function : Perl wrapper for generic commands module
## Returns  : @commands
## Arguments: $filehandle             => Filehandle to write to
##          : $infile_paths_ref       => Infile paths {REF}
##          : $orientation            => Orientation
##          : $outdir_path            => Out directory path
##          : $outfile_name           => Out file name
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdinfile_path         => Stdinfile path
##          : $stdoutfile_path        => Stdoutfile path
##          : $write_filenames        => Write filenames

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $infile_paths_ref;
    my $orientation;
    my $outdir_path;
    my $outfile_name;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdinfile_path;
    my $stdoutfile_path;
    my $write_filenames;

    my $tmpl = {
        filehandle => {
            store => \$filehandle,
        },
        infile_paths_ref => {
            required    => 1,
            default     => [],
            store       => \$infile_paths_ref,
            strict_type => 1,
        },
        orientation => {
            allow       => [qw{ landscape portrait }],
            required    => 1,
            store       => \$orientation,
            strict_type => 1,
        },
        outdir_path => {
            required    => 1,
            store       => \$outdir_path,
            strict_type => 1,
        },
        outfile_name => {
            required    => 1,
            store       => \$outfile_name,
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
        write_filenames => {
            allow       => [ undef, 0, 1 ],
            store       => \$write_filenames,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = ( get_executable_base_command( { base_command => $BASE_COMMAND, } ) );

    push @commands, join $SPACE, map { q{--infile} . $SPACE . $_ } @{$infile_paths_ref};

    push @commands, q{--orientation} . $SPACE . $orientation;

    push @commands, q{--outfolder} . $SPACE . $outdir_path;

    push @commands, q{--outfile} . $SPACE . $outfile_name;

    if ($write_filenames) {

        push @commands, q{--write-filenames};
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
