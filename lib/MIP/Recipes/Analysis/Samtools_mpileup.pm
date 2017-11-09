package MIP::Recipes::Analysis::Samtools_mpileup;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use autodie qw{ :all };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use File::Spec::Functions qw{ catdir catfile };

## CPANM
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_samtools_mpileup };

}

## Constants
Readonly my $ASTERISK   => q{*};
Readonly my $UNDERSCORE => q{_};
Readonly my $DOT        => q{.};
Readonly my $SPACE      => q{ };
Readonly my $PIPE       => q{|};
Readonly my $DASH       => q{-};
Readonly my $NEWLINE    => qq{\n};

sub analysis_samtools_mpileup {

## Function : samtools_mpileup
## Returns  :
## Arguments: $parameter_href          => Parameter hash {REF}
##          : $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $program_name            => Program name
##          : $family_id               => Family id
##          : $temp_directory          => Temporary directory
##          : $outaligner_dir          => Outaligner_dir used in the analysis
##          : $call_type               => Variant call type
##          : $xargs_file_counter      => The xargs file counter

  my ($arg_href) = @_;

  ## Default(s)
  my $family_id;
  my $outaligner_dir;
  my $temp_directory;
  my $call_type;
  my $xargs_file_counter;

  ## Flatten argument(s)
  my $parameter_href;
  my $active_parameter_href;
  my $sample_info_href;
  my $file_info_href;
  my $infile_lane_prefix_href;
  my $job_id_href;
  my $program_name;

  my $tmpl = {
      parameter_href => {
          required    => 1,
          defined     => 1,
          default     => {},
          strict_type => 1,
          store       => \$parameter_href,
      },
      active_parameter_href => {
          required    => 1,
          defined     => 1,
          default     => {},
          strict_type => 1,
          store       => \$active_parameter_href,
      },
      sample_info_href => {
          required    => 1,
          defined     => 1,
          default     => {},
          strict_type => 1,
          store       => \$sample_info_href,
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
      program_name => {
          required    => 1,
          defined     => 1,
          strict_type => 1,
          store       => \$program_name,
      },
      family_id => {
          default     => $arg_href->{active_parameter_href}{family_id},
          strict_type => 1,
          store       => \$family_id,
      },
      outaligner_dir => {
          default     => $arg_href->{active_parameter_href}{outaligner_dir},
          strict_type => 1,
          store       => \$outaligner_dir,
      },
      call_type =>
        { default => q{BOTH}, strict_type => 1, store => \$call_type, },
      xargs_file_counter => {
          default     => 0,
          allow       => qr/ ^\d+$ /xsm,
          strict_type => 1,
          store       => \$xargs_file_counter,
      },
  };

  check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  use MIP::File::Format::Pedigree qw{ create_fam_file };
  use MIP::Get::File qw{ get_file_suffix get_merged_infile_prefix };
  use MIP::IO::Files qw{ migrate_file xargs_migrate_contig_files };
  use MIP::Processmanagement::Slurm_processes
    qw{ slurm_submit_job_sample_id_dependency_add_to_family };
  use MIP::Program::Alignment::Samtools qw{ samtools_mpileup };
  use MIP::Program::Variantcalling::Bcftools
    qw{ bcftools_call bcftools_filter bcftools_norm };
  use MIP::Program::Variantcalling::Gatk qw{ gatk_concatenate_variants };
  use MIP::Program::Variantcalling::Perl qw{ replace_iupac };
  use MIP::QC::Record qw{ add_program_outfile_to_sample_info };
  use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
  use MIP::Script::Setup_script qw{ setup_script };
  use MIP::Set::File qw{ set_file_suffix };

  ## Retrieve logger object
  my $log = Log::Log4perl->get_logger(q{MIP});

  ## Set MIP program name
  my $mip_program_name = q{p} . $program_name;
  my $mip_program_mode = $active_parameter_href->{$mip_program_name};

  ## Unpack parameters
  my $job_id_chain = $parameter_href->{$mip_program_name}{chain};
  my $core_number =
    $active_parameter_href->{module_core_number}{$mip_program_name};
  my $time = $active_parameter_href->{module_time}{$mip_program_name};

  ## Alias
  my $program_outdirectory_name =
    $parameter_href->{ $mip_program_name }{outdir_name};
  my $program_directory = catfile( $outaligner_dir, $program_outdirectory_name );
  my $xargs_file_path_prefix;

  ## Filehandles
  # Create anonymous filehandle
  my $FILEHANDLE = IO::Handle->new();
  my $XARGSFILEHANDLE = IO::Handle->new();

  ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
  my ( $file_path, $program_info_path ) = setup_script(
      {
          active_parameter_href => $active_parameter_href,
          job_id_href           => $job_id_href,
          FILEHANDLE            => $FILEHANDLE,
          directory_id          => $family_id,
          program_name          => $program_name,
          program_directory => $program_directory,
          core_number => $core_number,
          process_time => $time,
          temp_directory => $temp_directory,
      }
  );

  ## Assign directories
  my $outfamily_file_directory =
    catfile( $active_parameter_href->{outdata_dir},
      $family_id, $outaligner_dir, $program_outdirectory_name );

  my $outfamily_directory = catfile( $active_parameter_href->{outdata_dir},
        $family_id, $outaligner_dir, $program_outdirectory_name );

  #Used downstream
  $parameter_href->{ $mip_program_name }{indirectory} = $outfamily_directory;

  ## Assign file_tags
  my $outfile_tag =
    $file_info_href->{$family_id}{ $mip_program_name }{file_tag};
  my $outfile_prefix = $family_id . $outfile_tag . $call_type;
  my $outfile_path_prefix = catfile( $temp_directory, $outfile_prefix );

  ## Assign suffix
  # Get infile_suffix from baserecalibration jobid chain
  my $infile_suffix = get_file_suffix(
      {
          parameter_href => $parameter_href,
          suffix_key     => q{alignment_file_suffix},
          jobid_chain    => $parameter_href->{pgatk_baserecalibration}{chain},
      }
  );

  ## Set file suffix for next module within jobid chain
  my $outfile_suffix = set_file_suffix(
      {
          parameter_href => $parameter_href,
          suffix_key     => q{variant_file_suffix},
          job_id_chain   => $job_id_chain,
          file_suffix => $parameter_href->{ $mip_program_name }{outfile_suffix},
      }
  );

  my $fam_file_path = catfile( $outfamily_file_directory, $family_id . q{.fam} );
  ## Create .fam file to be used in variant calling analyses
  create_fam_file(
      {
          parameter_href        => $parameter_href,
          active_parameter_href => $active_parameter_href,
          sample_info_href      => $sample_info_href,
          FILEHANDLE            => $FILEHANDLE,
          fam_file_path => $fam_file_path,
          include_header => 0,
      }
  );

  my %file_path_prefix;

SAMPLE:
  foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } )
  {
      ## Assign directories
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
      my $infile_tag =
        $file_info_href->{$sample_id}{pgatk_baserecalibration}{file_tag};
      my $infile_prefix = $merged_infile_prefix . $infile_tag;

      $file_path_prefix{$sample_id} =
        catfile( $temp_directory, $infile_prefix );

      ## Copy file(s) to temporary directory
      say {$FILEHANDLE} q{## Copy file(s) to temporary directory};

      ($xargs_file_counter) = xargs_migrate_contig_files(
          {
              FILEHANDLE      => $FILEHANDLE,
              XARGSFILEHANDLE => $XARGSFILEHANDLE,
              contigs_ref => \@{ $file_info_href->{contigs_size_ordered} },
              file_path   => $file_path,
              program_info_path  => $program_info_path,
              core_number        => $core_number,
              xargs_file_counter => $xargs_file_counter,
              infile             => $infile_prefix,
              indirectory        => $insample_directory,
              file_ending        => substr( $infile_suffix, 0, 2 )
                . $ASTERISK,
              temp_directory => $temp_directory,
          }
      );

      ## SamTools mpileup
      say {$FILEHANDLE} q{## SamTools mpileup};

      ## Create file commands for xargs
      ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
          {
              FILEHANDLE         => $FILEHANDLE,
              XARGSFILEHANDLE    => $XARGSFILEHANDLE,
              file_path          => $file_path,
              program_info_path  => $program_info_path,
              core_number        => $core_number,
              xargs_file_counter => $xargs_file_counter,
          }
      );

    }

      ## Split per contig
    CONTIG:
      foreach my $contig ( @{ $file_info_href->{contigs_size_ordered} } ) {

          ## Assemble file paths by adding file ending
          my @file_paths =
            map { $file_path_prefix{$_} . $UNDERSCORE . $contig . $infile_suffix }
            @{ $active_parameter_href->{sample_ids} };

          samtools_mpileup(
              {
                  infile_paths_ref => \@file_paths,
                  output_tags_ref  => [ qw{ DV AD } ],
                  referencefile_path =>
                    $active_parameter_href->{human_genome_reference},
                  stderrfile_path => $xargs_file_path_prefix . $DOT
                    . $contig
                    . q{.stderr.txt},
                  output_bcf                       => 1,
                  adjust_mq                        => 50,
                  per_sample_increased_sensitivity => 1,
                  FILEHANDLE                       => $XARGSFILEHANDLE,
              }
          );

          # Print pipe
          print {$XARGSFILEHANDLE} $PIPE . $SPACE;

          ## Get parameter
          my $samples_file;
          my $constrain;
          if ( $parameter_href->{dynamic_parameter}{trio} ) {

              $samples_file =
                catfile( $outfamily_file_directory, $family_id . q{.fam} )
                . $SPACE;
              $constrain = q{trio};
          }
          bcftools_call(
              {
                  form_fields_ref     => [ qw{ GQ }],
                  variants_only       => 1,
                  multiallelic_caller => 1,
                  samples_file        => $samples_file,
                  constrain           => $constrain,
                  stderrfile_path     => $xargs_file_path_prefix . $DOT
                    . $contig
                    . $UNDERSCORE
                    . q{call.stderr.txt},
                  FILEHANDLE => $XARGSFILEHANDLE,
              }
          );

          # Print pipe
          print {$XARGSFILEHANDLE} $PIPE . $SPACE;

          bcftools_filter(
              {
                  stderrfile_path => $xargs_file_path_prefix . $DOT
                    . $contig
                    . $UNDERSCORE
                    . q{filter.stderr.txt},
                  soft_filter => q{LowQual},
                  snp_gap     => 3,
                  indel_gap   => 10,
                  exclude =>
  q?\'%QUAL<10 || (RPB<0.1 && %QUAL<15) || (AC<2 && %QUAL<15) || %MAX(DV)<=3 || %MAX(DV)/%MAX(DP)<=0.25\'?,
                  FILEHANDLE => $XARGSFILEHANDLE,
              }
          );

          if ( $active_parameter_href->{replace_iupac} ) {

              ## Replace the IUPAC code in alternative allels with N for input stream and writes to stream
              replace_iupac(
                  {
                      stderrfile_path => $xargs_file_path_prefix . $DOT
                        . $contig
                        . $UNDERSCORE
                        . q{replace_iupac.stderr.txt},
                      FILEHANDLE => $XARGSFILEHANDLE
                  }
              );
          }

          # Print pipe
          print {$XARGSFILEHANDLE} $PIPE . $SPACE;

          ## BcfTools norm, Left-align and normalize indels, split multiallelics
          bcftools_norm(
              {
                  FILEHANDLE => $XARGSFILEHANDLE,
                  reference_path =>
                    $active_parameter_href->{human_genome_reference},
                  output_type  => q{v},
                  outfile_path => $outfile_path_prefix . $UNDERSCORE
                    . $contig
                    . $outfile_suffix,
                  multiallelic    => $DASH,
                  stderrfile_path => catfile(
                          $xargs_file_path_prefix . $DOT
                        . $contig
                        . $UNDERSCORE
                        . q{norm.stderr.txt}
                  ),
              }
          );
          say {$XARGSFILEHANDLE} $NEWLINE;

  }

  ## Writes sbatch code to supplied filehandle to concatenate variants in vcf format. Each array element is combined with the infile prefix and postfix.
  gatk_concatenate_variants(
      {
          active_parameter_href => $active_parameter_href,
          FILEHANDLE            => $FILEHANDLE,
          elements_ref          => \@{ $file_info_href->{contigs} },
          infile_prefix         => $outfile_path_prefix . $UNDERSCORE,
          infile_postfix        => $outfile_suffix,
          outfile_path_prefix   => $outfile_path_prefix,
          outfile_suffix        => $outfile_suffix,
      }
  );

  ## Copies file from temporary directory.
  say {$FILEHANDLE} q{## Copy file from temporary directory};
  migrate_file(
      {
          infile_path  => $outfile_path_prefix . $outfile_suffix . $ASTERISK,
          outfile_path => $outfamily_directory,
          FILEHANDLE   => $FILEHANDLE,
      }
  );
  say {$FILEHANDLE} q{wait}, $NEWLINE;

  close $FILEHANDLE;
  close $XARGSFILEHANDLE;

  if ( $active_parameter_href->{ $mip_program_name } == 1 ) {

      ## Collect samtools version in qccollect
      add_program_outfile_to_sample_info(
          {
              sample_info_href => $sample_info_href,
              program_name     => q{samtools},
              outdirectory     => $outfamily_directory,
              outfile          => $outfile_prefix . $outfile_suffix,
          }
      );
      ## Locating samtools_mpileup file
      add_program_outfile_to_sample_info(
          {
              sample_info_href => $sample_info_href,
              program_name     => q{samtools_mpileup},
              outdirectory     => $outfamily_directory,
              outfile          => $outfile_prefix . $outfile_suffix,
          }
      );
      add_program_outfile_to_sample_info(
          {
              sample_info_href => $sample_info_href,
              program_name     => q{bcftools},
              outdirectory     => $outfamily_directory,
              outfile          => $outfile_prefix . $outfile_suffix,
          }
      );

      slurm_submit_job_sample_id_dependency_add_to_family(
          {
              job_id_href             => $job_id_href,
              infile_lane_prefix_href => $infile_lane_prefix_href,
              sample_ids_ref   => \@{ $active_parameter_href->{sample_ids} },
              family_id        => $family_id,
              path             => $job_id_chain,
              log              => $log,
              sbatch_file_name => $file_path,
          }
      );
  }
  return;
}

1;
