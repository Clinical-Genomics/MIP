package MIP::Recipes::Analysis::Variant_integrity;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile devnull splitpath };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
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
    our $VERSION = 1.07;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ analysis_variant_integrity };

}

## Constants
Readonly my $ASTERISK   => q{*};
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

sub analysis_variant_integrity {

## Function : Tests sample for correct relatives (only performed for samples with relatives defined in pedigree file) performed on sequence data.
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $case_id               => Family id
##          : $file_info_href          => The file_info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $recipe_name            => Program name
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $temp_directory          => Temporary directory

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

    use MIP::File::Format::Pedigree qw{ create_fam_file };
    use MIP::Get::File qw{ get_io_files };
    use MIP::Get::Parameter qw{ get_recipe_parameters get_recipe_attributes };
    use MIP::Parse::File qw{ parse_io_outfiles };
    use MIP::Processmanagement::Processes qw{ submit_recipe };
    use MIP::Program::Variantcalling::Variant_integrity
      qw{ variant_integrity_mendel variant_integrity_father };
    use MIP::QC::Record qw{ add_recipe_metafile_to_sample_info };
    use MIP::Script::Setup_script qw{ setup_script };

    ### PREPROCESSING:

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Unpack parameters
    ## Get the io infiles per chain and id
    my %io = get_io_files(
        {
            id             => $case_id,
            file_info_href => $file_info_href,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => q{in},
            temp_directory => $temp_directory,
        }
    );
    my $infile_name_prefix = $io{in}{file_name_prefix};
    my $infile_path        = $io{in}{file_path};

    my $job_id_chain = get_recipe_attributes(
        {
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            attribute      => q{chain},
        }
    );
    my $recipe_mode = $active_parameter_href->{$recipe_name};
    my @sample_ids  = @{ $active_parameter_href->{sample_ids} };
    my $is_trio     = $parameter_href->{dynamic_parameter}{trio};
    my ( $core_number, $time, @source_environment_cmds ) = get_recipe_parameters(
        {
            active_parameter_href => $active_parameter_href,
            recipe_name           => $recipe_name,
        }
    );

    ## Set outfiles depending on sample data
    my %var_int_outfiles = (
        trio   => { mendel => q{txt}, },
        father => { father => q{txt}, },
    );

    ## Check if a father is included in analysis
    my $has_father;
  SAMPLE_ID:
    foreach my $sample_id (@sample_ids) {

        # Alias
        my $is_father =
          $sample_info_href->{sample}{$sample_id}{father};

        if ($is_father) {
            $has_father = 1;
        }
    }

    my @var_int_outfiles;
  MODE:
    while ( my ( $mode, $program_href ) = each %var_int_outfiles ) {

      VAR_INT_PROGRAM:
        while ( my ( $file_name_prefix, $file_suffix ) = each %{$program_href} ) {

            if ( $is_trio and $mode eq q{trio} ) {

                push @var_int_outfiles, $file_name_prefix . $DOT . $file_suffix;
                next;
            }
            if ( $has_father and $mode eq q{father} ) {

                push @var_int_outfiles, $file_name_prefix . $DOT . $file_suffix;
            }
        }
    }

    ## Check if we ahve something to analyse
    if ( not scalar @var_int_outfiles ) {

        $log->warn(q{Not a trio nor a father in analysis - skipping 'variant_integrity'});
        return;
    }
    ## Set and get the io files per chain, id and stream
    %io = (
        %io,
        parse_io_outfiles(
            {
                chain_id         => $job_id_chain,
                id               => $case_id,
                file_info_href   => $file_info_href,
                file_name_prefix => $infile_name_prefix,
                iterators_ref    => \@var_int_outfiles,
                outdata_dir      => $active_parameter_href->{outdata_dir},
                parameter_href   => $parameter_href,
                recipe_name      => $recipe_name,
                temp_directory   => $temp_directory,
            }
        )
    );

    my $outdir_path_prefix = $io{out}{dir_path_prefix};
    my %outfile_path       = %{ $io{out}{file_path_href} };

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    ## Creates recipe directories (info & data & script), recipe script filenames and writes sbatch header
    my ( $recipe_file_path, $recipe_info_path ) = setup_script(
        {
            active_parameter_href           => $active_parameter_href,
            core_number                     => $core_number,
            directory_id                    => $case_id,
            FILEHANDLE                      => $FILEHANDLE,
            job_id_href                     => $job_id_href,
            log                             => $log,
            process_time                    => $time,
            recipe_directory                => $recipe_name,
            recipe_name                     => $recipe_name,
            source_environment_commands_ref => \@source_environment_cmds,
        }
    );

    my $case_file_path = catfile( $outdir_path_prefix, $case_id . $DOT . q{fam} );

    ## Create .fam file to be used in variant calling analyses
    create_fam_file(
        {
            active_parameter_href => $active_parameter_href,
            fam_file_path         => $case_file_path,
            FILEHANDLE            => $FILEHANDLE,
            parameter_href        => $parameter_href,
            sample_info_href      => $sample_info_href,
        }
    );

    ## Variant_integrity
    ## If a trio
    if ($is_trio) {

        variant_integrity_mendel(
            {
                case_file    => $case_file_path,
                case_type    => $active_parameter_href->{genmod_models_case_type},
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => $infile_path,
                outfile_path => $outfile_path{mendel},
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    # If father is included in analysis
    if ($has_father) {

        variant_integrity_father(
            {
                case_file    => $case_file_path,
                case_type    => $active_parameter_href->{genmod_models_case_type},
                FILEHANDLE   => $FILEHANDLE,
                infile_path  => $infile_path,
                outfile_path => $outfile_path{father},
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }

    close $FILEHANDLE;

    if ( $recipe_mode == 1 ) {

        while ( my ( $outfile_tag, $outfile_path ) = each %outfile_path ) {

## Collect QC metadata info for later use
            add_recipe_metafile_to_sample_info(
                {
                    metafile_tag     => $outfile_tag,
                    path             => $outfile_path,
                    recipe_name      => $recipe_name,
                    sample_info_href => $sample_info_href,
                }
            );
        }

        submit_recipe(
            {
                dependency_method       => q{case_to_island},
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
    return;
}

1;
