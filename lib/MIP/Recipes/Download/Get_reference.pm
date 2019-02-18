package MIP::Recipes::Download::Get_reference;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
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
use MIP::Constants qw{ $NEWLINE $SPACE $UNDERSCORE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ get_reference };
}

sub get_reference {

## Function : Write get reference recipe (download, decompress and validate)
## Returns  :
## Arguments: $FILEHANDLE     => Filehandle to write to
##          : $parameter_href => Holds all parameters {REF}
##          : $recipe_name    => Recipe name
##          : $reference_href => Reference hash {REF}
##          : $quiet          => Quiet (no output)
##          : $verbose        => Verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $parameter_href;
    my $recipe_name;
    my $reference_href;

    ## Default(s)
    my $quiet;
    my $verbose;

    my $tmpl = {
        FILEHANDLE     => { defined => 1, required => 1, store => \$FILEHANDLE, },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
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
            default     => $arg_href->{parameter_href}{verbose},
            store       => \$verbose,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::File qw{ check_file_md5sum };
    use MIP::File::Decompression qw{ decompress_file };
    use MIP::Program::Download::Wget qw{ wget };

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
        my $outfile_path = catfile( $parameter_href->{reference_dir}, $outfile );
        my $url          = $reference_href->{url_prefix} . $file;

        ## Download
        say {$FILEHANDLE} q{## Download } . $recipe_name . $NEWLINE;

        wget(
            {
                FILEHANDLE   => $FILEHANDLE,
                outfile_path => $outfile_path,
                quiet        => $quiet,
                url          => $url,
                verbose      => $verbose,
            }
        );
        say {$FILEHANDLE} $NEWLINE;

        ## Check if file needs to be decompress and write decompression if so
        decompress_file(
            {
                FILEHANDLE   => $FILEHANDLE,
                outdir_path  => $parameter_href->{reference_dir},
                outfile_path => $outfile_path,
                decompress_program =>
                  $reference_href->{ q{out} . $key . $UNDERSCORE . q{decompress} },
            }
        );

        ## Check file integrity of file
        check_file_md5sum(
            {
                FILEHANDLE    => $FILEHANDLE,
                md5_file_path => $outfile_path,
                check_method =>
                  $reference_href->{ q{out} . $key . $UNDERSCORE . q{method} },
            }
        );
    }
    return;
}

1;
