package MIP::Recipes::Analysis::Picardtools_genotypeconcordance;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
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

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_picardtools_genotypeconcordance };

}

## Constants
Readonly my $ASTERISK   => q{*};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub analysis_picardtools_genotypeconcordance {

## Function : Compare metrics for this analysis run with the NIST reference dataset.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $family_id               => Family id
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $sample_id               => Sample id
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $reference_dir           => MIP reference directory
##          : $temp_directory          => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $program_name;
    my $sample_id;
    my $sample_info_href;

    ## Default(s)
    my $family_id;
    my $reference_dir;
    my $temp_directory;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        family_id => {
            default     => $arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name,
        },
        reference_dir => {
            default     => $arg_href->{active_parameter_href}{reference_dir},
            strict_type => 1,
            store       => \$reference_dir,
        },
        sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id,
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        temp_directory => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_module_parameters get_program_attributes };
    use MIP::Gnu::Coreutils qw{ gnu_cat };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_family_dead_end };
    use MIP::Program::Interval::Picardtools qw{ picardtools_intervallisttools };
    use MIP::Program::Variantcalling::Bcftools
      qw{ bcftools_norm bcftools_rename_vcf_samples bcftools_stats bcftools_rename_vcf_samples };
    use MIP::Program::Variantcalling::Gatk qw{ gatk_selectvariants };
    use MIP::Program::Variantcalling::Picardtools
      qw{ picardtools_genotypeconcordance };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Return if not nist_id
    return if ( not $sample_id =~ /$active_parameter_href->{nist_id}/sxm );

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Unpack parameters
    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            id             => $family_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            program_name   => $program_name,
            stream         => q{in},
            temp_directory => $temp_directory,
        }
    );
    my $infile_name_prefix = $io{in}{file_name_prefix};
    my $infile_path        = $io{in}{file_path};

    my $job_id_chain = get_program_attributes(
        {
            parameter_href => $parameter_href,
            program_name   => $program_name,
            attribute      => q{chain},
        }
    );
    my $nist_id            = $active_parameter_href->{nist_id};
    my $program_mode       = $active_parameter_href->{$program_name};
    my $referencefile_path = $active_parameter_href->{human_genome_reference};
    my ( $core_number, $time, @source_environment_cmds ) =
      get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            program_name          => $program_name,
        }
      );

    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id         => $job_id_chain,
                id               => $family_id,
                file_info_href   => $file_info_href,
                file_name_prefix => $infile_name_prefix,
                iterators_ref    => [$nist_id],
                outdata_dir      => $active_parameter_href->{outdata_dir},
                parameter_href   => $parameter_href,
                program_name     => $program_name,
                temp_directory   => $temp_directory,
            }
        )
    );

    my $outdir_path_prefix  = $io{out}{dir_path_prefix};
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my $outfile_suffix      = $io{out}{file_suffix};

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ($file_path) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $core_number,
            directory_id          => $family_id,
            error_trap => 0,    # Special case to allow "vcf.idx" to be created
            FILEHANDLE        => $FILEHANDLE,
            job_id_href       => $job_id_href,
            log               => $log,
            process_time      => $time,
            program_directory => $program_name,
            program_name      => $program_name,
            set_errexit => 0,    # Special case to allow "vcf.idx" to be created
            source_environment_commands_ref => \@source_environment_cmds,
            temp_directory                  => $temp_directory,
        }
    );

    my $nist_file_path = catfile( $temp_directory, q{NIST} );

    ## Rename vcf samples. The samples array will replace the sample names in the same order as supplied.
    bcftools_rename_vcf_samples(
        {
            FILEHANDLE => $FILEHANDLE,
            index      => 0,
            infile => $active_parameter_href->{nist_high_confidence_call_set},
            outfile_path_prefix => $nist_file_path . $UNDERSCORE . q{refrm},
            output_type         => q{v},
            temp_directory      => $temp_directory,
            sample_ids_ref      => [ $nist_id . q{-NIST} ],
        }
    );

    ## Modify since different ref genomes
    say {$FILEHANDLE} q{## Modify since different ref genomes};

    ## Do not print GL contigs
    print {$FILEHANDLE}
      q?perl -nae 'unless($_=~/##contig=<ID=GL\d+/) {print $_}' ?;

    ## Infile
    print {$FILEHANDLE} $nist_file_path . $UNDERSCORE . q{refrm.vcf} . $SPACE;

    ## Outfile
    print {$FILEHANDLE} q{>}
      . $SPACE
      . $nist_file_path
      . $DOT . q{vcf}
      . $SPACE;
    say {$FILEHANDLE} $NEWLINE;

    ## Bcftools stats
    say {$FILEHANDLE} q{## bcftools stats};
    bcftools_stats(
        {
            FILEHANDLE      => $FILEHANDLE,
            infile_path     => $nist_file_path . $DOT . q{vcf},
            stdoutfile_path => $nist_file_path . $DOT . q{vcf.stats},
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Generate ".idx" for downstream Picard by failling this process
    say {$FILEHANDLE}
      q{## Generate '.idx' for downstream Picard by failling this process};

    gatk_selectvariants(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $nist_file_path . $DOT . q{vcf},
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            memory_allocation  => q{Xmx2g},
            outfile_path       => $nist_file_path . q{XXX.vcf},
            referencefile_path => $referencefile_path,
            sample_names_ref   => [ $sample_id . q{XXX} ],
            temp_directory     => $temp_directory,
            verbosity          => $active_parameter_href->{gatk_logging_level},
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Create .interval_list file from NIST union bed
    say {$FILEHANDLE} q{## Prepare .interval_list file from NIST union bed},
      $NEWLINE;

    my $genome_dict_file_path = catfile( $temp_directory,
            $file_info_href->{human_genome_reference_name_prefix}
          . $DOT
          . q{dict} );
    ## Do not print contigs
    print {$FILEHANDLE}
q?perl -nae 'unless($_=~/NC_007605/ || $_=~/hs37d5/ || $_=~/GL\d+/) {print $_}' ?;
    print {$FILEHANDLE} catfile( $reference_dir,
        $file_info_href->{human_genome_reference_name_prefix} . $DOT . q{dict} )
      . $SPACE;
    print {$FILEHANDLE} q{>} . $SPACE . $genome_dict_file_path . $SPACE;
    say   {$FILEHANDLE} $NEWLINE;

    gnu_cat(
        {
            FILEHANDLE       => $FILEHANDLE,
            infile_paths_ref => [
                $genome_dict_file_path,
                $active_parameter_href->{nist_high_confidence_call_set_bed}
            ],
            outfile_path => $nist_file_path . $DOT . q{bed.dict_body},
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE}
q{## Remove target annotations, 'track', 'browse' and keep only 5 columns};
    print {$FILEHANDLE}
q?perl  -nae 'if ($_=~/@/) {print $_;} elsif ($_=~/^track/) {} elsif ($_=~/^browser/) {} else {print @F[0], "\t", (@F[1] + 1), "\t", @F[2], "\t", "+", "\t", "-", "\n";}' ?;
    ## Infile
    print {$FILEHANDLE} $nist_file_path . $DOT . q{bed.dict_body} . $SPACE;

    ## Remove unnecessary info and reformat
    print {$FILEHANDLE} q{>}
      . $SPACE
      . $nist_file_path
      . $DOT
      . q{bed.dict_body_col_5.interval_list}
      . $SPACE;
    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE} q{## Create }
      . $active_parameter_href->{nist_high_confidence_call_set_bed}
      . $DOT
      . q{interval_list};

    picardtools_intervallisttools(
        {
            FILEHANDLE => $FILEHANDLE,
            infile_paths_ref =>
              [ $nist_file_path . $DOT . q{bed.dict_body_col_5.interval_list} ],
            java_jar => catfile(
                $active_parameter_href->{picardtools_path},
                q{picard.jar}
            ),
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            memory_allocation  => q{Xmx2g},
            outfile_path       => $nist_file_path . $DOT . q{bed.interval_list},
            referencefile_path => $referencefile_path,
            temp_directory     => $active_parameter_href->{temp_directory},
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ### MIP data
    my $base_file_path = $outfile_path_prefix . $UNDERSCORE . $sample_id;
    ## GATK SelectVariants
    say {$FILEHANDLE} q{## GATK SelectVariants};

    my $select_outfile_path = $base_file_path . $outfile_suffix;
    gatk_selectvariants(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $infile_path,
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            memory_allocation  => q{Xmx2g},
            outfile_path       => $select_outfile_path,
            referencefile_path => $referencefile_path,
            sample_names_ref   => [$sample_id],
            temp_directory     => $temp_directory,
            verbosity          => $active_parameter_href->{gatk_logging_level},
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Left align, normalize and split allels
    say {$FILEHANDLE} q{## bcftools norm};

    my $norm_outfile_path =
      $base_file_path . $UNDERSCORE . q{norm} . $outfile_suffix;
    bcftools_norm(
        {
            FILEHANDLE     => $FILEHANDLE,
            infile_path    => $select_outfile_path,
            multiallelic   => q{-},
            outfile_path   => $norm_outfile_path,
            reference_path => $referencefile_path,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Modify since different ref genomes
    say {$FILEHANDLE} q{## Modify since different ref genomes};
    my $reformat_outfile_path =
      $base_file_path . $UNDERSCORE . q{norm_refrm.vcf};

    print {$FILEHANDLE}
q?perl -nae 'unless($_=~/##contig=<ID=NC_007605,length=171823>/ || $_=~/##contig=<ID=hs37d5,length=35477943>/ || $_=~/##contig=<ID=GL\d+/) {print $_}' ?;

    ## Infile
    print {$FILEHANDLE} $norm_outfile_path . $SPACE;

    ## Outfile
    print {$FILEHANDLE} q{>} . $SPACE . $reformat_outfile_path . $SPACE;
    say   {$FILEHANDLE} $NEWLINE;

    ## BcfTools Stats
    say {$FILEHANDLE} q{## bcftools stats};
    bcftools_stats(
        {
            FILEHANDLE      => $FILEHANDLE,
            infile_path     => $reformat_outfile_path,
            stdoutfile_path => $reformat_outfile_path . $DOT . q{stats},
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Generate ".idx" for downstream Picard by failling this process
    say {$FILEHANDLE}
      q{## Generate '.idx' for downstream Picard by failling this process};

    gatk_selectvariants(
        {
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $reformat_outfile_path,
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            memory_allocation  => q{Xmx2g},
            outfile_path       => $reformat_outfile_path . q{XXX.vcf},
            referencefile_path => $referencefile_path,
            sample_names_ref   => [ $sample_id . q{XXX} ],
            temp_directory     => $temp_directory,
            verbosity          => $active_parameter_href->{gatk_logging_level},
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE}
q{## Picard GenotypeConcordance - Genome restricted by union - good quality};
    picardtools_genotypeconcordance(
        {
            call_sample   => $sample_id,
            FILEHANDLE    => $FILEHANDLE,
            infile_path   => $reformat_outfile_path,
            intervals_ref => [ $nist_file_path . $DOT . q{bed.interval_list} ],
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            java_jar => catfile(
                $active_parameter_href->{picardtools_path},
                q{picard.jar}
            ),
            memory_allocation    => q{Xmx2g},
            min_depth            => 10,
            min_genotype_quality => 20,
            outfile_prefix_path  => $base_file_path . $UNDERSCORE . q{bed},
            referencefile_path   => $referencefile_path,
            temp_directory       => $active_parameter_href->{temp_directory},
            truth_file_path      => $nist_file_path . $DOT . q{vcf},
            truth_sample         => $nist_id . q{-NIST},
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    say {$FILEHANDLE} q{## Picard GenotypeConcordance - Genome - good quality};
    picardtools_genotypeconcordance(
        {
            call_sample => $sample_id,
            FILEHANDLE  => $FILEHANDLE,
            infile_path => $reformat_outfile_path,
            java_use_large_pages =>
              $active_parameter_href->{java_use_large_pages},
            java_jar => catfile(
                $active_parameter_href->{picardtools_path},
                q{picard.jar}
            ),
            memory_allocation    => q{Xmx2g},
            min_depth            => 10,
            min_genotype_quality => 20,
            outfile_prefix_path  => $base_file_path,
            referencefile_path   => $referencefile_path,
            truth_file_path      => $nist_file_path . $DOT . q{vcf},
            truth_sample         => $nist_id . q{-NIST},
            temp_directory       => $active_parameter_href->{temp_directory},
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    close $FILEHANDLE;

    if ( $program_mode == 1 ) {

        slurm_submit_job_sample_id_dependency_family_dead_end(
            {
                family_id               => $family_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                path                    => $job_id_chain,
                sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
                sbatch_file_name => $file_path,
            }
        );
    }
    return;
}

1;
