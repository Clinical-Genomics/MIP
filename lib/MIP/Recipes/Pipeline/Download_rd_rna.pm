package MIP::Recipes::Pipeline::Download_rd_rna;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
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
    our $VERSION = 1.04;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ pipeline_download_rd_rna };
}

sub pipeline_download_rd_rna {

## Function : Download references recipes for rd_rna pipeline
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this download hash {REF}
##          : $quiet                 => Be quiet
##          : $temp_directory        => Temporary directory
##          : $verbose               => Verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;

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
        quiet => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
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

    use MIP::Recipes::Download::1000g_indels qw{ download_1000g_indels };
    use MIP::Recipes::Download::Ctat_resource_lib qw{ download_ctat_resource_lib };
    use MIP::Recipes::Download::Dbsnp qw{ download_dbsnp };
    use MIP::Recipes::Download::Gencode_annotation qw{ download_gencode_annotation };
    use MIP::Recipes::Download::Get_reference qw{ get_reference };
    use MIP::Recipes::Download::Human_reference qw{ download_human_reference };
    use MIP::Recipes::Download::Mills_and_1000g_indels
      qw{ download_mills_and_1000g_indels };
    use MIP::Recipes::Download::Pfam qw{ download_pfam };
    use MIP::Recipes::Download::Runstatus qw{ download_runstatus };

    ## Retrieve logger object now that log_file has been set
    my $log = Log::Log4perl->get_logger( uc q{mip_download} );

    ### Download recipes
    ## Create code reference table for download recipes
    my %download_recipe = (
        q{1000g_indels}        => \&download_1000g_indels,
        ctat_resource_lib      => \&download_ctat_resource_lib,
        dbsnp                  => \&download_dbsnp,
        gencode_annotation     => \&download_gencode_annotation,
        human_reference        => \&download_human_reference,
        mills_and_1000g_indels => \&download_mills_and_1000g_indels,
        pfam                   => \&download_pfam,
    );

    # Storing job_ids from SLURM. However, all recipes are currently indepent
    my %job_id;

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

                ## Save paths for run status check downstream
                push @{ $active_parameter_href->{runstatus_paths} }, $outfile_path;

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
                        genome_version        => $genome_version,
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

    return if ( not @{ $active_parameter_href->{runstatus_paths} } );

    ## Check downloaded file exists and has a file size greater than zero
    download_runstatus(
        {
            active_parameter_href => $active_parameter_href,
            job_id_href           => \%job_id,
            recipe_name           => q{runstatus},
            temp_directory        => $temp_directory,
        }
    );

    return;
}

1;
