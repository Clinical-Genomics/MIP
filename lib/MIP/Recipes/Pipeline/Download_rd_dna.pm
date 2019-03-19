package MIP::Recipes::Pipeline::Download_rd_dna;

use 5.026;
use Carp;
use charnames qw{ :full :short };
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
    our @EXPORT_OK = qw{ pipeline_download_rd_dna };
}

sub pipeline_download_rd_dna {

## Function : Download references recipes for rd_dna pipeline
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

    use MIP::Recipes::Download::1000g_all_sv qw{ download_1000g_all_sv };
    use MIP::Recipes::Download::1000g_all_wgs qw{ download_1000g_all_wgs };
    use MIP::Recipes::Download::1000g_indels qw{ download_1000g_indels };
    use MIP::Recipes::Download::1000g_omni qw{ download_1000g_omni };
    use MIP::Recipes::Download::1000g_sites qw{ download_1000g_sites };
    use MIP::Recipes::Download::1000g_snps qw{ download_1000g_snps };
    use MIP::Recipes::Download::Clinvar qw{ download_clinvar };
    use MIP::Recipes::Download::Dbnsfp qw{ download_dbnsfp };
    use MIP::Recipes::Download::Dbsnp qw{ download_dbsnp };
    use MIP::Recipes::Download::Genomic_superdups qw{ download_genomic_superdups };
    use MIP::Recipes::Download::Get_reference qw{ get_reference };
    use MIP::Recipes::Download::Giab qw{ download_giab };
    use MIP::Recipes::Download::Gnomad qw{ download_gnomad };
    use MIP::Recipes::Download::Hapmap qw{ download_hapmap };
    use MIP::Recipes::Download::Human_reference qw{ download_human_reference };
    use MIP::Recipes::Download::Mills_and_1000g_indels
      qw{ download_mills_and_1000g_indels };

    ## Retrieve logger object now that log_file has been set
    my $log = Log::Log4perl->get_logger( uc q{mip_download} );

    ### Download recipes
    ## Create code reference table for download recipes
    my %download_recipe = (
        q{1000g_all_sv}        => \&download_1000g_all_sv,
        q{1000g_all_wgs}       => \&download_1000g_all_wgs,
        q{1000g_indels}        => \&download_1000g_indels,
        q{1000g_omni}          => \&download_1000g_omni,
        q{1000g_sites}         => \&download_1000g_sites,
        q{1000g_snps}          => \&download_1000g_snps,
        clinvar                => \&download_clinvar,
        dbnsfp                 => \&download_dbnsfp,
        dbsnp                  => \&download_dbsnp,
        genomic_superdups      => \&download_genomic_superdups,
        giab                   => \&download_giab,
        gnomad                 => \&download_gnomad,
        hapmap                 => \&download_hapmap,
        human_reference        => \&download_human_reference,
        mills_and_1000g_indels => \&download_mills_and_1000g_indels,
    );

    # Storing job_ids from SLURM, however currently all are independent
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

                ## Check if reference file already exists in reference directory
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
                        genome_version        => $genome_version,
                        job_id_href           => \%job_id,
                        quiet                 => $quiet,
                        recipe_name           => $reference_id,
                        reference_href        => $reference_href,
                        reference_version     => $reference_version,
                        temp_directory        => $temp_directory,
                    }
                );
            }
        }
    }

    return;
}

1;
