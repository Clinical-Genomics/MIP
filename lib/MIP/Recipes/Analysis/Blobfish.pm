package MIP::Recipes::Analysis::Blobfish;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile splitpath };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use List::MoreUtils qw { uniq };
use Readonly;

## MIPs lib/
use MIP::Constants
  qw{ $ASTERISK $AMPERSAND $COLON $DOT $LOG_NAME $NEWLINE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.07;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_blobfish };

}

sub analysis_blobfish {

## Function : Blobfish (DESeq2) recipe for wts data
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File info hash {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
    my $profile_base_command;
    my $temp_directory;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        case_id => {
            default     => $arg_href->{active_parameter_href}{case_id},
            store       => \$case_id,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
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
        profile_base_command => {
            default     => q{sbatch},
            store       => \$profile_base_command,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Blobfish qw{ blobfish_allvsall };
    use MIP::Sample_info qw{ set_file_path_to_store set_recipe_outfile_in_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Collect infiles for all sample_ids
    my @sample_indir_paths;
    my @sample_phenotypes;
  SAMPLE_ID:
    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        ## Get the io infiles per chain and id
        my %sample_io = get_io_files(
            {
                id             => $sample_id,
                file_info_href => $file_info_href,
                parameter_href => $parameter_href,
                recipe_name    => $recipe_name,
                stream         => q{in},
                temp_directory => $temp_directory,
            }
        );
        push @sample_indir_paths, $sample_io{in}{dir_path};

        ## Get condition
        push @sample_phenotypes, $sample_info_href->{sample}{$sample_id}{phenotype};
    }

    ## Get outdir_path
    my $outdir_path =
      catdir( $active_parameter_href->{outdata_dir}, $case_id, $recipe_name );

    ## Get recipe attributes and parameters
    my $job_id_chain = get_recipe_attributes(
        {
            attribute      => q{chain},
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
        }
    );
    my $recipe_mode     = $active_parameter_href->{$recipe_name};
    my %recipe_resource = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Filehandles
    # Create anonymous filehandle
    my $filehandle = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $recipe_resource{core_number},
            directory_id                    => $case_id,
            filehandle                      => $filehandle,
            job_id_href                     => $job_id_href,
            memory_allocation               => $recipe_resource{memory},
            process_time                    => $recipe_resource{time},
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            temp_directory                  => $temp_directory,
        }
    );

    ### SHELL

    say {$filehandle} q{## Generate tx2gene file};
    my $gtf_file_path     = $active_parameter_href->{transcript_annotation};
    my $tx2gene_file_path = catfile( $outdir_path, q{tx2gene.txt} );
    _generate_tx2gene_file(
        {
            filehandle        => $filehandle,
            gtf_file_path     => $gtf_file_path,
            tx2gene_file_path => $tx2gene_file_path,
        }
    );

    say {$filehandle} q{## BlobFish};
    blobfish_allvsall(
        {
            conditions_ref    => \@sample_phenotypes,
            filehandle        => $filehandle,
            indir_paths_ref   => \@sample_indir_paths,
            outdir_path       => $outdir_path,
            tx2gene_file_path => $tx2gene_file_path,
        }
    );

    close $filehandle;

    if ( $recipe_mode == 1 ) {

        ## Collect QC metadata info for later use
        my @conditions = uniq @sample_phenotypes;
        my $de_outfile_name =
            $conditions[0]
          . $UNDERSCORE . q{vs}
          . $UNDERSCORE
          . $conditions[1]
          . q{.results.tsv};
        set_recipe_outfile_in_sample_info(
            {
                path             => catfile( $outdir_path, $de_outfile_name ),
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

        set_file_path_to_store(
            {
                format           => q{meta},
                id               => $case_id,
                path             => catfile( $outdir_path, $de_outfile_name ),
                recipe_name      => $recipe_name,
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                base_command         => $profile_base_command,
                dependency_method    => q{case_to_island},
                case_id              => $case_id,
                job_id_chain         => $job_id_chain,
                job_id_href          => $job_id_href,
                job_reservation_name => $active_parameter_href->{job_reservation_name},
                log                  => $log,
                max_parallel_processes_count_href =>
                  $file_info_href->{max_parallel_processes_count},
                recipe_file_path   => $recipe_file_path,
                sample_ids_ref     => \@{ $active_parameter_href->{sample_ids} },
                submission_profile => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

sub _generate_tx2gene_file {

## Function : Generate tx2gene file for Blobfish
## Returns  :
## Arguments: $filehandle       => Filehandle
##          : $gtf_file_path    => Path to annotation file
##          : tx2gene_file_path => Outfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $gtf_file_path;
    my $tx2gene_file_path;

    my $tmpl = {
        filehandle => {
            required => 1,
            store    => \$filehandle,
        },
        gtf_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$gtf_file_path,
            strict_type => 1,
        },
        tx2gene_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$tx2gene_file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

     use MIP::Environment::Executable qw{ get_executable_base_command };

    my @commands = ( get_executable_base_command( { base_command => q{perl}, } ), );

    # Execute perl
    my $tx2gene_generator = join $SPACE, @commands;
    ## Print header and initiate hash
    $tx2gene_generator .=
      q? -nae 'BEGIN {print q{TXNAME,GENEID} . qq{\n}; %txgene;}?;
    ## When the file has been processed; print hash
    $tx2gene_generator .=
      q? END {foreach $tx (keys %txgene){print $tx . q{,} . $txgene{$tx} .qq{\n}; } }?;
    ## Skip headers
    $tx2gene_generator .= q? if ($_ =~/^#/){next;}?;
    ## Check for keywords in attribute field
    $tx2gene_generator .= q? if (($F[8] =~ /gene_id/) and ($F[10] =~ /transcript_id/))?;
    ## Capture gene, transcript id and remove [;"] from the names
    $tx2gene_generator .=
      q? {$gene = $F[9]; $gene =~ tr/[;"]//d; $tx = $F[11]; $tx =~tr/[;"]//d;?;
    ## Store in hash
    $tx2gene_generator .= q? $txgene{$tx} = $gene;} else{next;}'?;

    say {$filehandle} $tx2gene_generator
      . $SPACE
      . $gtf_file_path
      . $SPACE . q{>}
      . $SPACE
      . $tx2gene_file_path;

    return;
}

1;
