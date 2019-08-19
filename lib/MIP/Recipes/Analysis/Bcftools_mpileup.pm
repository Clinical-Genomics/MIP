package MIP::Recipes::Analysis::Bcftools_mpileup;

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
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $ASTERISK $DASH $DOT $NEWLINE $PIPE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.11;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_bcftools_mpileup };

}

## Constants
Readonly my $ADJUST_MQ => 50;
Readonly my $SNP_GAP   => 3;
Readonly my $INDEL_GAP => 10;

sub analysis_bcftools_mpileup {

## Function : Call snv and indels using Bcftools mpileup
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id                 => Family id
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $profile_base_command    => Submission profile base command
##          : $recipe_name             => Program name
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory
##          : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_href;
    my $recipe_name;
    my $sample_info_href;

    ## Default(s)
    my $case_id;
    my $profile_base_command;
    my $temp_directory;
    my $xargs_file_counter;

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
        temp_directory_ref => {
            default     => $arg_href->{active_parameter_href}{temp_directory},
            store       => \$temp_directory,
            strict_type => 1,
        },
        xargs_file_counter => {
            allow       => qr/ ^\d+$ /xsm,
            default     => 0,
            store       => \$xargs_file_counter,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Format::Pedigree qw{ create_fam_file };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_attributes get_recipe_resources };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Variantcalling::Bcftools
      qw{ bcftools_call bcftools_filter bcftools_mpileup bcftools_norm bcftools_view };
    use MIP::Program::Variantcalling::Gatk qw{ gatk_concatenate_variants };
    use MIP::Program::Variantcalling::Perl qw{ replace_iupac };
    use MIP::Sample_info qw{ set_recipe_outfile_in_sample_info };
    use MIP::Recipes::Analysis::Xargs qw{ xargs_command };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger( uc q{mip_analyse} );

    ## Unpack parameters
    my $job_id_chain = get_recipe_attributes(
        {
            attribute      => q{chain},
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
        }
    );
    my $recipe_mode     = $active_parameter_href->{$recipe_name};
    my $reference_path  = $active_parameter_href->{human_genome_reference};
    my %recipe_resource = get_recipe_resources(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );
    my $core_number = $recipe_resource{core_number};

    my $xargs_file_path_prefix;

    ## Set and get the io files per chain, id and stream
    my %io = parse_io_outfiles(
        {
            chain_id               => $job_id_chain,
            id                     => $case_id,
            file_info_href         => $file_info_href,
            file_name_prefixes_ref => [$case_id],
            outdata_dir            => $active_parameter_href->{outdata_dir},
            parameter_href         => $parameter_href,
            recipe_name            => $recipe_name,
            temp_directory         => $temp_directory,
        }
    );
    my $outdir_path_prefix  = $io{out}{dir_path_prefix};
    my $outfile_path_prefix = $io{out}{file_path_prefix};
    my $outfile_suffix      = $io{out}{file_suffix};
    my $outfile_path        = $outfile_path_prefix . $outfile_suffix;

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE      = IO::Handle->new();
    my $XARGSFILEHANDLE = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $case_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            log                             => $log,
            memory_allocation               => $recipe_resource{memory},
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            process_time                    => $recipe_resource{time},
            source_environment_commands_ref => $recipe_resource{load_env_ref},
            temp_directory                  => $temp_directory,
        }
    );

    ## Create .fam file to be used in variant calling analyses
    my $fam_file_path = catfile( $outdir_path_prefix, $case_id . $DOT . q{fam} );
    create_fam_file(
        {
            active_parameter_href => $active_parameter_href,
            fam_file_path         => $fam_file_path,
            FILEHANDLE            => $FILEHANDLE,
            include_header        => 0,
            log                   => $log,
            parameter_href        => $parameter_href,
            sample_info_href      => $sample_info_href,
        }
    );

    ### SHELL:

    ## Collect infiles for all sample_ids to enable migration to temporary directory
    my %mpileup_infile_path;
    while ( my ( $sample_id_index, $sample_id ) =
        each @{ $active_parameter_href->{sample_ids} } )
    {

        ## Get the io infiles per chain and id
        my %sample_io = get_io_files(
            {
                id             => $sample_id,
                file_info_href => $file_info_href,
                parameter_href => $parameter_href,
                recipe_name    => $recipe_name,
                stream         => q{in},
            }
        );
        my $infile_suffix      = $sample_io{in}{file_suffix};
        my $infile_path_prefix = $sample_io{in}{file_path_prefix};
        my $infile_path        = $infile_path_prefix . $infile_suffix;

        ## Store temp infile path for each sample_id and contig
        $mpileup_infile_path{$sample_id} =
          $sample_io{in}{file_path_href};
    }

    ## Bcftools mpileup
    say {$FILEHANDLE} q{## Bcftools mpileup};

    ## Create file commands for xargs
    ( $xargs_file_counter, $xargs_file_path_prefix ) = xargs_command(
        {
            core_number        => $core_number,
            FILEHANDLE         => $FILEHANDLE,
            file_path          => $recipe_file_path,
            recipe_info_path   => $recipe_info_path,
            XARGSFILEHANDLE    => $XARGSFILEHANDLE,
            xargs_file_counter => $xargs_file_counter,
        }
    );

  CONTIG:
    foreach my $contig ( @{ $file_info_href->{contigs_size_ordered} } ) {

        ## Assemble contig file paths for mpileup
        my @mpileup_file_paths =
          map { $mpileup_infile_path{$_}{$contig} }
          @{ $active_parameter_href->{sample_ids} };

        my $stderrfile_path_prefix = $xargs_file_path_prefix . $DOT . $contig;
        my $output_type            = q{b};

        bcftools_mpileup(
            {
                adjust_mq                        => $ADJUST_MQ,
                FILEHANDLE                       => $XARGSFILEHANDLE,
                infile_paths_ref                 => \@mpileup_file_paths,
                output_tags_ref                  => [qw{ DV AD }],
                output_type                      => $output_type,
                per_sample_increased_sensitivity => 1,
                referencefile_path               => $reference_path,
                regions_ref                      => [$contig],
                stderrfile_path => $stderrfile_path_prefix . $DOT . q{stderr.txt},
            }
        );

        # Print pipe
        print {$XARGSFILEHANDLE} $PIPE . $SPACE;

        ## Get parameter
        my $samples_file;
        ## Special case: bcftools version 1.9 does not output GQ when constraint is applied
        my $constrain = $active_parameter_href->{bcftools_mpileup_constrain};
        if ( $parameter_href->{cache}{trio} and $constrain ) {

            $samples_file =
              catfile( $outdir_path_prefix, $case_id . $DOT . q{fam} ) . $SPACE;
            $constrain = q{trio};
        }

        bcftools_call(
            {
                constrain           => $constrain,
                FILEHANDLE          => $XARGSFILEHANDLE,
                form_fields_ref     => [qw{ GQ }],
                multiallelic_caller => 1,
                output_type         => $output_type,
                samples_file_path   => $samples_file,
                stderrfile_path     => $stderrfile_path_prefix
                  . $UNDERSCORE
                  . q{call.stderr.txt},
                variants_only => 1,
            }
        );

        if ( $active_parameter_href->{bcftools_mpileup_filter_variant} ) {

            # Print pipe
            print {$XARGSFILEHANDLE} $PIPE . $SPACE;

            bcftools_filter(
                {
                    FILEHANDLE      => $XARGSFILEHANDLE,
                    exclude         => _build_bcftools_filter_expr(),
                    indel_gap       => $INDEL_GAP,
                    output_type     => $output_type,
                    snp_gap         => $SNP_GAP,
                    soft_filter     => q{LowQual},
                    stderrfile_path => $stderrfile_path_prefix
                      . $UNDERSCORE
                      . q{filter.stderr.txt},
                }
            );
        }

        if ( not $active_parameter_href->{bcftools_mpileup_keep_unnormalised} ) {

            # Print pipe
            print {$XARGSFILEHANDLE} $PIPE . $SPACE;

            bcftools_norm(
                {
                    FILEHANDLE      => $XARGSFILEHANDLE,
                    multiallelic    => $DASH,
                    output_type     => $output_type,
                    reference_path  => $reference_path,
                    stderrfile_path => $stderrfile_path_prefix
                      . $UNDERSCORE
                      . q{norm.stderr.txt},
                }
            );
        }

        # Print pipe
        print {$XARGSFILEHANDLE} $PIPE . $SPACE;

        my $bcftools_outfile_path;

        if ( not $active_parameter_href->{replace_iupac} ) {

            ## End stream and write to disc
            $bcftools_outfile_path =
              $outfile_path_prefix . $DOT . $contig . $outfile_suffix;
        }
        bcftools_view(
            {
                FILEHANDLE      => $XARGSFILEHANDLE,
                outfile_path    => $bcftools_outfile_path,
                output_type     => q{v},
                stderrfile_path => $stderrfile_path_prefix
                  . $UNDERSCORE
                  . q{view.stderr.txt},
            }
        );
        if ( $active_parameter_href->{replace_iupac} ) {

            print {$XARGSFILEHANDLE} $PIPE . $SPACE;

            ## End stream and write to disc
            $bcftools_outfile_path =
              $outfile_path_prefix . $DOT . $contig . $outfile_suffix;
            ## Replace the IUPAC code in alternative allels with N for input stream and writes to stream
            replace_iupac(
                {
                    FILEHANDLE      => $XARGSFILEHANDLE,
                    stderrfile_path => $stderrfile_path_prefix
                      . $UNDERSCORE
                      . q{replace_iupac.stderr.txt},
                    stdoutfile_path => $bcftools_outfile_path,
                }
            );
        }

        say {$XARGSFILEHANDLE} $NEWLINE;
    }

    ## Writes sbatch code to supplied filehandle to concatenate variants in vcf format. Each array element is combined with the infile prefix and postfix.
    my $concat_infile_path_prefix = $outfile_path_prefix . $UNDERSCORE . q{ordered};
    gatk_concatenate_variants(
        {
            active_parameter_href => $active_parameter_href,
            FILEHANDLE            => $FILEHANDLE,
            elements_ref          => \@{ $file_info_href->{contigs} },
            infile_postfix        => $outfile_suffix,
            infile_prefix         => $outfile_path_prefix,
            outfile_path_prefix   => $outfile_path_prefix,
            outfile_suffix        => $outfile_suffix,
        }
    );

    close $FILEHANDLE or $log->logcroak(q{Could not close FILEHANDLE});
    close $XARGSFILEHANDLE
      or $log->logcroak(q{Could not close XARGSFILEHANDLE});

    if ( $recipe_mode == 1 ) {

        ## Locating bcftools_mpileup file
        set_recipe_outfile_in_sample_info(
            {
                path             => $outfile_path,
                recipe_name      => q{bcftools_mpileup},
                sample_info_href => $sample_info_href,
            }
        );

        submit_recipe(
            {
                base_command            => $profile_base_command,
                dependency_method       => q{sample_to_case},
                case_id                 => $case_id,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                log                     => $log,
                job_id_chain            => $job_id_chain,
                recipe_file_path        => $recipe_file_path,
                sample_ids_ref          => \@{ $active_parameter_href->{sample_ids} },
                submission_profile      => $active_parameter_href->{submission_profile},
            }
        );
    }
    return 1;
}

sub _build_bcftools_filter_expr {

## Function : Create filter expression for bcftools
## Returns  : $regexp

    ## Constants
    Readonly my $FILTER_SEPARATOR => q{ || };

    # Add minimum value for QUAL field
    my $expr = q?\'\"%QUAL<10?;

    # Add read position bias threshold
    $expr .= $FILTER_SEPARATOR . q{(RPB<0.1 && %QUAL<15)};

    # Add allele count expression
    $expr .= $FILTER_SEPARATOR . q{(AC<2 && %QUAL<15)};

    # Add number of high-qual non-reference bases
    $expr .= $FILTER_SEPARATOR . q{%MAX(DV)<=3};

    # Add high-qual non-reference bases / high-qual bases
    $expr .= $FILTER_SEPARATOR . q?%MAX(DV)/%MAX(DP)<=0.25\"\'?;

    return $expr;
}

1;
