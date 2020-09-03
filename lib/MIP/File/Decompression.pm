package MIP::File::Decompression;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

## MIPs lib/
use MIP::Constants qw{ $DOT $NEWLINE $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.05;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ decompress_files };
}

sub decompress_files {

## Function : Decompress file(s) with gzip, unzip or tar
## Returns  : @commands
##          : $file_path      => Path to file
##          : $file_paths_ref => path to file
##          : $filehandle     => Filehandle to write to
##          : $outdir_path    => Outdirectory path
##          : $outfile_path   => Outfile path
##          : $program        => Decompress the file(s) using program supplied

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_path;
    my $file_paths_ref;
    my $filehandle;
    my $outdir_path;
    my $outfile_path;
    my $program;

    my $tmpl = {
        file_path => {
            defined     => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        file_paths_ref => {
            default     => $arg_href->{file_paths_ref} ||= [ $arg_href->{file_path} ],
            defined     => 1,
            required    => 1,
            store       => \$file_paths_ref,
            strict_type => 1,
        },
        filehandle => {
            defined => 1,
            store   => \$filehandle,
        },
        outdir_path => {
            default     => cwd(),
            store       => \$outdir_path,
            strict_type => 1,
        },
        outfile_path => {
            store       => \$outfile_path,
            strict_type => 1,
        },
        program => {
            allow       => [ undef, qw{ gzip tar unzip } ],
            store       => \$program,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Path qw{ remove_file_path_suffix };
    use MIP::Program::Gzip qw{ gzip };
    use MIP::Program::Tar qw{ tar };
    use MIP::Program::Zip qw{ unzip };

    return if ( not defined $program );

    croak q{Supply infile for tar!} if ( $program eq q{tar} and not $file_path );

    my %decompress_api = (
        gzip => {
            method   => \&gzip,
            arg_href => {
                decompress       => 1,
                filehandle       => $filehandle,
                force            => 1,
                infile_paths_ref => $file_paths_ref,
                outfile_path     => remove_file_path_suffix(
                    {
                        file_path         => $outfile_path,
                        file_suffixes_ref => [qw{ .gz }],
                    }
                ),
                quiet  => 1,
                stdout => 1,
            },
        },
        unzip => {
            method   => \&unzip,
            arg_href => {
                filehandle       => $filehandle,
                force            => 1,
                infile_paths_ref => $file_paths_ref,
                outdir_path      => $outdir_path,
            },
        },
        tar => {
            method   => \&tar,
            arg_href => {
                extract           => 1,
                file_path         => $file_path,
                filehandle        => $filehandle,
                filter_gzip       => 1,
                outdirectory_path => $outdir_path,
            },
        },
    );

    my @commands =
      $decompress_api{$program}{method}->( { %{ $decompress_api{$program}{arg_href} } } );

    if ( defined $filehandle ) {
        say {$filehandle} $NEWLINE;
    }

    return @commands;
}

1;
