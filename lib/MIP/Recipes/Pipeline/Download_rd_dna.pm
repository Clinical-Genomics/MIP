package MIP::Recipes::Pipeline::Download_rd_dna;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use List::MoreUtils qw { any };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $NEWLINE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ pipeline_download_rd_dna };
}

sub pipeline_download_rd_dna {

## Function : Download references recipes for rd_dna pipeline
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this download hash {REF}
##          : $FILEHANDLE            => Filehandle to write to
##          : $quiet                 => Be quiet
##          : $temp_directory        => Temporary directory
##          : $verbose               => Verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $FILEHANDLE;

    ## Default(s)
    my $quiet;
    my $temp_directory;
    my $verbose;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        FILEHANDLE => { defined => 1, required => 1, store => \$FILEHANDLE, },
        quiet      => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$quiet,
            strict_type => 1,
        },
        temp_directory => {
            defined     => 1,
            required    => 1,
            store       => \$temp_directory,
            strict_type => 1,
        },
        verbose => {
            default     => $arg_href->{active_parameter_href}{verbose},
            store       => \$verbose,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Gnu::Bash qw{ gnu_cd };
    use MIP::Gnu::Coreutils qw{ gnu_mkdir };
    use MIP::Recipes::Download::Get_reference qw{ get_reference };
    use MIP::Recipes::Download::Human_reference qw{ download_human_reference };

    ## Retrieve logger object now that log_file has been set
    my $log = Log::Log4perl->get_logger( uc q{mip_download} );

    my $pwd = cwd();

    ### Download recipes
    ## Create code reference table for download recipes
    my %download_recipe = (
        human_reference  => \&download_human_reference,
        q{1000g_all_wgs} => \&download_1000g_all_wgs,
    );

    # Storing job_ids from SLURM
    my %job_id;

    say {$FILEHANDLE} q{## Create reference directory};
    gnu_mkdir(
        {
            indirectory_path => $active_parameter_href->{reference_dir},
            parents          => 1,
            FILEHANDLE       => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Since all commands should assume working directory to be the reference directory
    gnu_cd(
        {
            directory_path => $active_parameter_href->{reference_dir},
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

  REFERENCE:
    while ( my ( $reference_id, $versions_ref ) =
        each %{ $active_parameter_href->{reference} } )
    {

        next REFERENCE if ( not exists $download_recipe{$reference_id} );

      REFERENCE_VERSION:
        foreach my $reference_version ( @{$versions_ref} ) {

          GENOME_VERSION:
            foreach my $genome_version (
                @{ $active_parameter_href->{reference_genome_versions} } )
            {

                my $reference_href =
                  $active_parameter_href->{reference_feature}{$reference_id}
                  {$genome_version}{$reference_version};

                next GENOME_VERSION
                  if (
                    not exists $active_parameter_href->{reference_feature}{$reference_id}
                    {$genome_version} );

                next GENOME_VERSION
                  if (
                    not exists $active_parameter_href->{reference_feature}{$reference_id}
                    {$genome_version}{$reference_version} );

                ## Build file name and path
                my $outfile_name = $reference_href->{outfile};
                my $outfile_path =
                  catfile( $active_parameter_href->{reference_dir}, $outfile_name );

                ## Check if reference already exists in reference directory
                next GENOME_VERSION if ( -f $outfile_path );

                $log->info( q{Cannot find reference file:} . $outfile_path );
                $log->info(
                        q{Will try to download: }
                      . $reference_id
                      . q{ version: }
                      . $reference_version,
                );

                $download_recipe{$reference_id}->(
                    {
                        active_parameter_href => $active_parameter_href,
                        job_id_href           => \%job_id,
                        recipe_name           => $reference_id,
                        reference_href        => $reference_href,
                        reference_version     => $reference_version,
                        quiet                 => $quiet,
                        temp_directory        => $temp_directory,
                    }
                );
            }
        }
    }

    ## Move back to original
    gnu_cd(
        {
            directory_path => $pwd,
            FILEHANDLE     => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;
    return;
}

1;
