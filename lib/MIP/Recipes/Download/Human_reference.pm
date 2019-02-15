package MIP::Recipes::Download::Human_reference;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile devnull };
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
use MIP::Constants qw{ $COLON $NEWLINE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ download_human_reference };

}

sub download_human_reference {

## Function : Download human genome reference builds
## Returns  :
## Arguments: $job_id_href    => The job_id hash {REF}
##          : $parameter_href => Parameter hash {REF}
##          : $reference_href => Reference hash {REF}
##          : $recipe_name    => Recipe name
##          : $quiet          => Quiet (no output)
##          : $temp_directory => Temporary directory for recipe
##          : $verbose        => Verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $reference_href;

    ## Default(s)
my $quiet;
my $temp_directory;
my $verbose;

    my $tmpl = {
job_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$job_id_href,
            strict_type => 1,
},
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
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$reference_href,
},
quiet      => {
            default     => 1,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$quiet
},
temp_directory => {
            store       => \$temp_directory,
            strict_type => 1,
},
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Cwd;
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(uc q{mip_download});

    ## Unpack parameters
    my @reference_genome_versions = @{$parameter_href->{reference_genome_versions}};
    my $file_path = catfile(cwd(), q{test.sh});

    ## Filehandle(s)
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

#    open $FILEHANDLE, q{>}, $file_path
#  or croak q{Cannot write to} . $SPACE . $file_path . $COLON . $SPACE . $OS_ERROR;

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $parameter_href,
            directory_id                    => q{mip_download},
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            log                             => $log,
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
        }
);

    ### SHELL:

    say {$FILEHANDLE} q{## } . $recipe_name;

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

 #       download(
 #           {
 #               parameter_href => $parameter_href,
 #               FILEHANDLE     => $FILEHANDLE,
 #               url            => $reference_href->{url_prefix} . $file,
 #               outfile_path   => $outfile_path,
 #               file_id        => $recipe_name,
 #           }
 #       );

        ## Check if file needs to be decompress and write decompression if so
#        decompress_file(
#            {
#                parameter_href => $parameter_href,
#                FILEHANDLE     => $FILEHANDLE,
#                outfile_path   => $outfile_path,
#                file_decompress =>
#                  $reference_href->{ q{out} . $key . $UNDERSCORE . q{decompress#} },
 #           }
 #       );

        ## Check file integrity of file
#        check_file(
#            {
#                FILEHANDLE         => $FILEHANDLE,
#                outfile_path       => $outfile_path,
#                outfile_path_check => $outfile_path,
#                check_method =>
#                  $reference_href->{ q{out} . $key . $UNDERSCORE . q{method} },
#            }
#        );
      }

    ## Close FILEHANDLES
    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});

    return 1;
}

1;
