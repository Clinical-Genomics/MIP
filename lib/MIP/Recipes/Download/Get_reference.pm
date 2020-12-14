package MIP::Recipes::Download::Get_reference;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $NEWLINE $SPACE $UNDERSCORE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ get_reference };
}

## Constants
Readonly my $READ_TIMEOUT_SEC => 20;
Readonly my $TIMEOUT_SEC      => 20;
Readonly my $DOWNLOAD_TRIES   => 12;
Readonly my $WAIT_RETRY_SEC   => 300;

sub get_reference {

## Function : Write get reference recipe (download, decompress and validate)
## Returns  :
## Arguments: $filehandle     => Filehandle to write to
##          : $recipe_name    => Recipe name
##          : $reference_dir  => Reference directory
##          : $reference_href => Reference hash {REF}
##          : $quiet          => Quiet (no output)
##          : $verbose        => Verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $recipe_name;
    my $reference_dir;
    my $reference_href;

    ## Default(s)
    my $outdir_path;
    my $quiet;
    my $verbose;

    my $tmpl = {
        filehandle  => { defined => 1, required => 1, store => \$filehandle, },
        outdir_path => {
            default     => $arg_href->{reference_dir},
            store       => \$outdir_path,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        reference_dir => {
            defined     => 1,
            required    => 1,
            store       => \$reference_dir,
            strict_type => 1,
        },
        reference_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$reference_href,
            strict_type => 1,
        },
        quiet => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$quiet,
            strict_type => 1,
        },
        verbose => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$verbose,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Validate::File qw{ check_file_md5sum };
    use MIP::File::Decompression qw{ decompress_files };
    use MIP::Program::Wget qw{ wget };

    ## Potential download files
    my @file_keys = qw{ file file_check
      file_index file_index_check };

  REFERENCE_FILE:
    foreach my $key (@file_keys) {

        next REFERENCE_FILE
          if ( not exists $reference_href->{$key} );

        ## Install reference
        my $file         = $reference_href->{$key};
        my $outfile      = $reference_href->{ q{out} . $key };
        my $outfile_path = catfile( $reference_dir, $outfile );
        my $url          = $reference_href->{url_prefix} . $file;
        my $user         = $reference_href->{user};

        ## Download
        say {$filehandle} q{## Download } . $recipe_name . $NEWLINE;

        wget(
            {
                continue          => 1,
                filehandle        => $filehandle,
                outfile_path      => $outfile_path,
                read_timeout      => $READ_TIMEOUT_SEC,
                retry_connrefused => 1,
                timeout           => $TIMEOUT_SEC,
                tries             => $DOWNLOAD_TRIES,
                wait_retry        => $WAIT_RETRY_SEC,
                quiet             => $quiet,
                url               => $url,
                user              => $user,
                verbose           => $verbose,
            }
        );
        say {$filehandle} $NEWLINE;

        ## Check if file needs to be decompress and write decompression if so
        decompress_files(
            {
                filehandle  => $filehandle,
                outdir_path => $outdir_path,
                file_path   => $outfile_path,
                program =>
                  $reference_href->{ q{out} . $key . $UNDERSCORE . q{decompress} },
                outfile_path => $outfile_path,
            }
        );

        ## Check file integrity of file
        check_file_md5sum(
            {
                filehandle    => $filehandle,
                md5_file_path => $outfile_path,
                check_method =>
                  $reference_href->{ q{out} . $key . $UNDERSCORE . q{method} },
            }
        );
    }
    return 1;
}

1;
