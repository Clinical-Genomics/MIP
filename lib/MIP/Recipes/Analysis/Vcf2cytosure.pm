package MIP::Recipes::Analysis::Vcf2cytosure;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
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
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_vcf2cytosure };

}

## Constants
Readonly my $AMPERSAND  => q{&};
Readonly my $ASTERISK   => q{*};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub analysis_vcf2cytosure {

## Function : Convert VCF with structural variations to the “.CGH” format used by the CytoSure Interpret Software
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $bin_size                => Bin size
##          : $family_id               => Family id
##          : $file_info_href          => File_info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $outfamily_directory     => Out family directory
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
    my $outfamily_directory;
    my $parameter_href;
    my $program_name;
    my $sample_info_href;

    ## Default(s)
    my $bin_size;
    my $family_id;
    my $outaligner_dir;
    my $temp_directory;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        bin_size => {
            default     => $arg_href->{active_parameter_href}{tiddit_bin_size},
            strict_type => 1,
            store       => \$bin_size
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
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            store       => \$outaligner_dir,
            strict_type => 1,
        },
        outfamily_directory => {
            defined     => 1,
            required    => 1,
            store       => \$outfamily_directory,
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

    use MIP::Cluster qw{ get_core_number };
    use MIP::IO::Files qw{ migrate_file };
    use MIP::Get::File qw{ get_file_suffix get_merged_infile_prefix };
    use MIP::Get::Parameter qw{ get_module_parameters get_program_parameters };
    use MIP::Program::Variantcalling::Vcf2cytosure qw{ vcf2cytosure_convert };
    use MIP::Processmanagement::Processes qw{ print_wait };
    use MIP::Processmanagement::Slurm_processes
      qw{ slurm_submit_job_sample_id_dependency_add_to_family };
    use MIP::Program::Variantcalling::Bcftools qw{ bcftools_view };
    use MIP::Program::Variantcalling::Tiddit qw{ tiddit_coverage };
    use MIP::QC::Record
      qw{ add_program_outfile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };
    use MIP::Set::File qw{ set_file_suffix };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Set MIP program name
    my $mip_program_name = q{p} . $program_name;
    my $mip_program_mode = $active_parameter_href->{$mip_program_name};

    ## Unpack parameters
    my $job_id_chain = $parameter_href->{$mip_program_name}{chain};
    my $max_cores_per_node = $active_parameter_href->{max_cores_per_node};
    my $modifier_core_number =
      scalar( @{ $active_parameter_href->{sample_ids} } );
    my $program_outdirectory_name =
      $parameter_href->{$mip_program_name}{outdir_name};
    my ( $core_number, $time, $source_environment_cmd ) = get_module_parameters(
        {
            active_parameter_href => $active_parameter_href,
            mip_program_name      => $mip_program_name,
        }
    );

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    $core_number = get_core_number(
        {
            max_cores_per_node   => $max_cores_per_node,
            modifier_core_number => $modifier_core_number,
            module_core_number   => $core_number,
        }
    );

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $family_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            process_time                    => $time,
            program_directory               => catfile($outaligner_dir, $program_outdirectory_name),
            program_name                    => $program_name,
            source_environment_commands_ref => [$source_environment_cmd],
            temp_directory                  => $temp_directory,
        }
    );

    ## Assign file_tags
    my %file_path_prefix;
    my $infile_prefix;
    my $infile_tag;
    my $sample_outfile_prefix;
    my $outfile_tag;

    my $process_batches_count = 1;

    # Loop over all samples
    while ( my ( $sample_id_index, $sample_id ) =
        each @{ $active_parameter_href->{sample_ids} } )
    {
        say {$FILEHANDLE} $NEWLINE . q{######################################################};
        say {$FILEHANDLE} q{Sample id index:} . $sample_id_index . q{, sample id:} . $sample_id;

        say {$FILEHANDLE} q{Creating coverage file with tiddit -cov for sample} . $SPACE . $sample_id;

        # Using tiddit coverage, create coverage file from .bam file of this sample
        my $insample_directory = catdir( $active_parameter_href->{outdata_dir},
            $sample_id, $outaligner_dir );

        ## Add merged infile name prefix after merging all BAM files per sample_id
        my $merged_infile_prefix = get_merged_infile_prefix(
            {
                file_info_href => $file_info_href,
                sample_id      => $sample_id,
            }
        );

        ## Assign file_tags
        $infile_tag = $file_info_href->{$sample_id}{pgatk_baserecalibration}{file_tag};
        $infile_prefix = $merged_infile_prefix . $infile_tag;
        $outfile_tag = $file_info_href->{$family_id}{ptiddit}{file_tag};
        $sample_outfile_prefix = $merged_infile_prefix . $outfile_tag;

        ## Assign suffix
        my $infile_suffix = get_file_suffix(
            {
                jobid_chain    => $parameter_href->{pgatk_baserecalibration}{chain},
                parameter_href => $parameter_href,
                suffix_key     => q{alignment_file_suffix},
            }
        );

        ## Set file suffix for next module within jobid chain
        my $outfile_suffix = set_file_suffix(
            {
                file_suffix => $parameter_href->{ptiddit}{coverage_file_suffix},
                job_id_chain   => $job_id_chain,
                parameter_href => $parameter_href,
                suffix_key     => q{variant_file_suffix},
            }
        );

        #q{.bam} -> ".b*" for getting index as well
        my $infile_path = catfile( $insample_directory,
            $infile_prefix . substr( $infile_suffix, 0, 2 ) . $ASTERISK );

        $file_path_prefix{$sample_id}{in} =
          catfile( $temp_directory, $infile_prefix );
        $file_path_prefix{$sample_id}{out} =
          catfile( $temp_directory, $sample_outfile_prefix );

        $process_batches_count = print_wait(
            {
                FILEHANDLE            => $FILEHANDLE,
                max_process_number    => $core_number,
                process_batches_count => $process_batches_count,
                process_counter       => $sample_id_index,
            }
        );


        ## Copy file(s) to temporary directory
        say {$FILEHANDLE} q{## Copy bam file to temporary directory};
        migrate_file(
            {
                  FILEHANDLE   => $FILEHANDLE,
                  infile_path  => $infile_path,
                  outfile_path => $temp_directory,
            }
        );
        say {$FILEHANDLE} q{wait}, $NEWLINE;

        ## Tiddit coverage
        tiddit_coverage(
            {
                bin_size    => $bin_size,
                FILEHANDLE  => $FILEHANDLE,
                infile_path => $file_path_prefix{$sample_id}{in}
                  . $infile_suffix,
                outfile_path_prefix => $file_path_prefix{$sample_id}{out},
            }
        );
        say {$FILEHANDLE} $AMPERSAND . $SPACE . $NEWLINE;
        say {$FILEHANDLE} q{wait}, $NEWLINE;


        # Extract SV from this sample from merged SV VCF file
        say {$FILEHANDLE} q{Using bcftools_view to extract SV for sample} . $SPACE . $sample_id;

        $infile_tag = $file_info_href->{$sample_id}{psv_combinevariantcallsets}{file_tag};

        say {$FILEHANDLE} $NEWLINE;
              say {$FILEHANDLE} q{insample_directory:} . $insample_directory;
              say {$FILEHANDLE} q{merged_infile_prefix:} . $merged_infile_prefix;
              say {$FILEHANDLE} q{infile_tag:} . $infile_tag;
              say {$FILEHANDLE} q{infile_prefix:} . $infile_prefix;
              say {$FILEHANDLE} q{sample_outfile_prefix:} . $sample_outfile_prefix;
              say {$FILEHANDLE} q{infile_path:} . $infile_path;

        #$infile_path = catfile( $insample_directory, $infile_prefix ),
        #say {$FILEHANDLE} q{insample_directory:} . $insample_directory;
        #say {$FILEHANDLE} q{$infile_prefix:} .$infile_prefix;
        #say {$FILEHANDLE} q{VCF directory:} . $infile_path ;
        #say {$FILEHANDLE} q{infile tag from psv_combinevariantcallsets:} . $infile_tag;


        #bcftools_view(
        #    {
        #        FILEHANDLE      => $FILEHANDLE,
        #        infile_path     => $file_path_prefix . $infile_suffix,
        #        samples_ref     => [$sample_id],
        #        stdoutfile_path => $sample_exonic_file_path,
        #    }
        #);
        say {$FILEHANDLE} $NEWLINE;




































    }

    say {$FILEHANDLE} q{wait}, $NEWLINE;











































    #my %file_path_prefix;
    #my $outfile_tag =
    #  $file_info_href->{$family_id}{$mip_program_name}{file_tag};
    #my $outfile_prefix = $family_id . $outfile_tag ;
    #my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix );

    ## Assign suffix
    #my $infile_suffix = get_file_suffix(
    #    {
    #        jobid_chain    => $parameter_href->{psv_combinevariantcallsets}{chain},
    #        parameter_href => $parameter_href,
    #        suffix_key     => q{variant_file_suffix},
    #    }
    #);

    #my $process_batches_count = 1;

    ## Collect infiles for all sample_ids to enable migration to temporary directory
    #  while ( my ( $sample_id_index, $sample_id ) =
    #      each @{ $active_parameter_href->{sample_ids} } )
    #  {

          ## Assign directories
    #      my $insample_directory = catdir( $active_parameter_href->{outdata_dir},
    #          $sample_id, $outaligner_dir );

          ## Add merged infile name prefix after merging all BAM files per sample_id
    #      my $merged_infile_prefix = get_merged_infile_prefix(
    #          {
    #              file_info_href => $file_info_href,
    #              sample_id      => $sample_id,
    #          }
    #      );

          ## Assign file_tags
    #      my $infile_tag =
    #        $file_info_href->{$sample_id}{psv_combinevariantcallsets}{file_tag};
    #      my $infile_prefix         = $merged_infile_prefix . $infile_tag;
    #      my $sample_outfile_prefix = $merged_infile_prefix . $outfile_tag;


    #      my $infile_path = catfile( $insample_directory, $infile_prefix);

    #      $file_path_prefix{$sample_id}{in} =
    #        catfile( $temp_directory, $infile_prefix );
    #      $file_path_prefix{$sample_id}{out} =
    #        catfile( $temp_directory, $sample_outfile_prefix );

    #      say {$FILEHANDLE} $NEWLINE;
    #      say {$FILEHANDLE} q{insample_directory:} . $insample_directory;
    #      say {$FILEHANDLE} q{merged_infile_prefix:} . $merged_infile_prefix;
    #      say {$FILEHANDLE} q{infile_tag:} . $infile_tag;
    #      say {$FILEHANDLE} q{infile_prefix:} . $infile_prefix;
    #      say {$FILEHANDLE} q{sample_outfile_prefix:} . $sample_outfile_prefix;
    #      say {$FILEHANDLE} q{infile_path:} . $infile_path;

    #      $process_batches_count = print_wait(
    #         {
    #              FILEHANDLE            => $FILEHANDLE,
    #              max_process_number    => $core_number,
    #              process_batches_count => $process_batches_count,
    #              process_counter       => $sample_id_index,
    #          }
    #      );

          ## Copy file(s) to temporary directory
    #      say {$FILEHANDLE} q{## Copy file(s) to temporary directory};
    #      migrate_file(
    #          {
    #              FILEHANDLE   => $FILEHANDLE,
    #              infile_path  => $infile_path,
    #              outfile_path => $temp_directory,
    #          }
    #      );
    #  }
    #  say {$FILEHANDLE} q{wait}, $NEWLINE;

      # Restart counter
      #$process_batches_count = 1;

      ## Collect infiles for all sample_ids
      #while ( my ( $sample_id_index, $sample_id ) =
      #    each @{ $active_parameter_href->{sample_ids} } )
      #{

      #    $process_batches_count = print_wait(
      #        {
      #            FILEHANDLE            => $FILEHANDLE,
      #            max_process_number    => $core_number,
      #            process_batches_count => $process_batches_count,
      #            process_counter       => $sample_id_index,
      #        }
      #    );

          ## Vcf2cytosure convert
      #    vcf2cytosure_convert(
      #        {
      #            coverage_file => ,
      #            frequency   =>
      #            frequency_tag =>
      #            FILEHANDLE  => $FILEHANDLE,
      #            infile_path => $file_path_prefix{$sample_id}{in}
      #              . $infile_suffix,
      #            no_filter =>
      #            stdoutfile_path =>
      #            variant_size =>
      #            vcf_infile_path =>
      #        }
      #    );
      #    say {$FILEHANDLE} $AMPERSAND . $SPACE . $NEWLINE;
      #}
      #say {$FILEHANDLE} q{wait}, $NEWLINE;

      ## Copies file from temporary directory.
      #say {$FILEHANDLE} q{## Copy file from temporary directory};
      #migrate_file(
      #    {
      #        FILEHANDLE   => $FILEHANDLE,
      #        infile_path  => $outfile_path_prefix . $outfile_suffix . $ASTERISK,
      #        outfile_path => $outfamily_directory,
      #    }
      #);
      #say {$FILEHANDLE} q{wait}, $NEWLINE;

      #close $FILEHANDLE;

      if ( $mip_program_mode == 1 ) {

      #add_program_outfile_to_sample_info(
      #    {
      #        path => catfile(
      #            $outfamily_directory, $outfile_prefix . $outfile_suffix
      #        ),
      #        program_name     => q{vcf2cytosure},
      #        sample_info_href => $sample_info_href,
      #    }
      #);

      slurm_submit_job_sample_id_dependency_add_to_family(
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
