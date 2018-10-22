package MIP::Recipes::Analysis::Blobfish;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile splitpath };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use List::MoreUtils qw { uniq };
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_blobfish };

}

## Constants
Readonly my $ASTERISK   => q{*};
Readonly my $AMPERSAND  => q{&};
Readonly my $COLON      => q{:};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub analysis_blobfish {

## Function : Blobfish (DESeq2) recipe for wts data
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id               => Family id
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $program_name;
    my $sample_info_href;

    ## Default(s)
    my $family_id;
    my $temp_directory;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            store       => \$family_id,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        infile_lane_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_lane_prefix_href,
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
        program_name => {
            defined     => 1,
            required    => 1,
            store       => \$program_name,
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
    use MIP::Get::Parameter qw{ get_module_parameters get_program_attributes };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_family_dead_end };
    use MIP::Program::Variantcalling::Blobfish qw{ blobfish_allvsall };
    use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Collect infiles for all sample_ids
    my @sample_indir_paths;
    my @sample_phenotypes;
    foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

        ## Get the io infiles per chain and id
        my %sample_io = get_io_files(
            {
                id             => $sample_id,
                file_info_href => $file_info_href,
                parameter_href => $parameter_href,
                program_name   => $program_name,
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
      catdir( $active_parameter_href->{outdata_dir}, $family_id, $program_name );

    ## Get program attributes and parameters
    my $job_id_chain = get_program_attributes(
        {
            attribute      => q{chain},
            parameter_href => $parameter_href,
            program_name   => $program_name,
        }
    );
    my $program_mode = $active_parameter_href->{$program_name};
    my ( $core_number, $time, @source_environment_cmds ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
    );

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $family_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            log                             => $log,
            process_time                    => $time,
            program_directory               => $program_name,
            program_name                    => $program_name,
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );

    ### SHELL

    say {$FILEHANDLE} q{## Generate tx2gene file};
    ## Print header and initiate hash
    my $tx2gene_generator =
      q?perl -nae 'BEGIN {print q{TXNAME,GENEID} . qq{\n}; %txgene;}?;
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

    my $gtf_file = $active_parameter_href->{transcripts_file};

    my $tx2gene_file_path = catfile( $outdir_path, q{tx2gene.txt} );

    say {$FILEHANDLE} $tx2gene_generator
      . $SPACE
      . $gtf_file
      . $SPACE . q{>}
      . $SPACE
      . $tx2gene_file_path;

    say {$FILEHANDLE} q{## BlobFish};
    blobfish_allvsall(
        {
            conditions_ref    => \@sample_phenotypes,
            FILEHANDLE        => $FILEHANDLE,
            indir_paths_ref   => \@sample_indir_paths,
            outdir_path       => $outdir_path,
            tx2gene_file_path => $tx2gene_file_path,
        }
    );

    close $FILEHANDLE;

    if ( $program_mode == 1 ) {

        ## Collect QC metadata info for later use
        add_program_outfile_to_sample_info(
            {
                path             => $outdir_path,
                program_name     => $program_name,
                sample_info_href => $sample_info_href,
            }
        );

        slurm_submit_job_sample_id_dependency_family_dead_end(
            {
                family_id               => $family_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                path                    => $job_id_chain,
                sample_ids_ref          => \@{ $active_parameter_href->{sample_ids} },
                sbatch_file_name        => $file_path,
            }
        );
    }
    return;
}

1;
