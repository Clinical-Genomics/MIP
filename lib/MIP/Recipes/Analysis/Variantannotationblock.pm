package MIP::Recipes::Analysis::Variantannotationblock;

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
    our @EXPORT_OK = qw{ analysis_variantannotationblock };

}

## Constants
Readonly my $DOT        => q{.};
Readonly my $NEWLINE    => qq{\n};
Readonly my $TAB        => qq{\t};
Readonly my $UNDERSCORE => q{_};

sub analysis_variantannotationblock {

## Function : Run consecutive module
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $call_type               => Variant call type
##          : $family_id               => Family id
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $outaligner_dir          => Outaligner dir used in the analysis
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $xargs_file_counter      => The xargs file counter

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
    my $call_type;
    my $family_id;
    my $outaligner_dir;
    my $xargs_file_counter;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        call_type =>
          { default => q{BOTH}, strict_type => 1, store => \$call_type, },
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
        outaligner_dir => {
            default     => $arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir,
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
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$xargs_file_counter,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Recipes::Analysis::Endvariantannotationblock
      qw{ analysis_endvariantannotationblock_rio };
    use MIP::Recipes::Analysis::Frequency_filter
      qw{ analysis_frequency_filter_rio };
    use MIP::Recipes::Analysis::Mip_vcfparser qw{ analysis_mip_vcfparser_rio };
    use MIP::Recipes::Analysis::Prepareforvariantannotationblock
      qw{ analysis_prepareforvariantannotationblock_rio };
    use MIP::Recipes::Analysis::Rankvariant
      qw{ analysis_rankvariant_rio analysis_rankvariant_rio_unaffected };
    use MIP::Recipes::Analysis::Rhocall qw{ analysis_rhocall_annotate_rio };
    use MIP::Recipes::Analysis::Snpeff qw{ analysis_snpeff_rio };
    use MIP::Recipes::Analysis::Vep qw{ analysis_vep_rio };
    use MIP::Recipes::Analysis::Vt qw{ analysis_vt_rio };
    use MIP::Script::Setup_script qw{ setup_script };

    ## Constants
    Readonly my $PROCESS_TIME => 80;

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my $core_number = $active_parameter_href->{max_cores_per_node};

    ## Filehandles
    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    # Set order of supplying user info
    my @rio_program_order =
      qw{ prepareforvariantannotationblock rhocall vt frequency_filter varianteffectpredictor vcfparser snpeff rankvariant endvariantannotationblock };

    # Store what to supply to user
    my %rio_program = (
        prepareforvariantannotationblock =>
          q{[Prepareforvariantannotationblock]},
        rhocall                   => q{[rhocall]},
        vt                        => q{[Vt]},
        frequency_filter          => q{[Frequency filter]},
        varianteffectpredictor    => q{[Varianteffectpredictor]},
        vcfparser                 => q{[Vcfparser]},
        snpeff                    => q{[Snpeff]},
        rankvariant               => q{[Rankvariant]},
        endvariantannotationblock => q{[Endvariantannotationblock]},
    );

  RIO_PROGRAM:
    foreach my $program_name (@rio_program_order) {

        if ( $active_parameter_href->{$program_name} ) {

            my $program_header = $rio_program{$program_name};

            $log->info( $TAB . $program_header );
        }
    }

    ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
    my ( $file_path, $program_info_path ) = setup_script(
        {
            active_parameter_href => $active_parameter_href,
            core_number           => $core_number,
            directory_id          => $family_id,
            FILEHANDLE            => $FILEHANDLE,
            job_id_href           => $job_id_href,
            log                   => $log,
            process_time          => $PROCESS_TIME,
            program_directory     => $outaligner_dir,
            program_name          => $program_name,
        }
    );

    ## Copy files for variantannotationblock to enable restart and skip of modules within block
    if ( $active_parameter_href->{prepareforvariantannotationblock} ) {

        ($xargs_file_counter) = analysis_prepareforvariantannotationblock_rio(
            {
                active_parameter_href   => $active_parameter_href,
                call_type               => $call_type,
                FILEHANDLE              => $FILEHANDLE,
                file_info_href          => $file_info_href,
                file_path               => $file_path,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                parameter_href          => $parameter_href,
                program_info_path       => $program_info_path,
                program_name            => q{prepareforvariantannotationblock},
                sample_info_href        => $sample_info_href,
                stderr_path        => $program_info_path . $DOT . q{stderr.txt},
                xargs_file_counter => $xargs_file_counter,
            }
        );
    }
    if ( $active_parameter_href->{rhocall} ) {

        my $infamily_directory = catdir( $active_parameter_href->{outdata_dir},
            $family_id, $outaligner_dir );
        my $outfamily_directory = $infamily_directory;

        ($xargs_file_counter) = analysis_rhocall_annotate_rio(
            {
                active_parameter_href   => $active_parameter_href,
                call_type               => $call_type,
                FILEHANDLE              => $FILEHANDLE,
                file_info_href          => $file_info_href,
                file_path               => $file_path,
                infamily_directory      => $infamily_directory,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                outfamily_directory     => $outfamily_directory,
                parameter_href          => $parameter_href,
                program_info_path       => $program_info_path,
                program_name            => q{rhocall},
                sample_info_href        => $sample_info_href,
                stderr_path        => $program_info_path . $DOT . q{stderr.txt},
                xargs_file_counter => $xargs_file_counter,
            }
        );
    }
    if ( $active_parameter_href->{vt} ) {

        my $infamily_directory = catdir( $active_parameter_href->{outdata_dir},
            $family_id, $outaligner_dir );
        my $outfamily_directory = $infamily_directory;

        ($xargs_file_counter) = analysis_vt_rio(
            {
                active_parameter_href   => $active_parameter_href,
                call_type               => $call_type,
                FILEHANDLE              => $FILEHANDLE,
                file_info_href          => $file_info_href,
                file_path               => $file_path,
                infamily_directory      => $infamily_directory,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                parameter_href          => $parameter_href,
                outfamily_directory     => $outfamily_directory,
                program_info_path       => $program_info_path,
                program_name            => q{vt},
                sample_info_href        => $sample_info_href,
                stderr_path        => $program_info_path . $DOT . q{stderr.txt},
                xargs_file_counter => $xargs_file_counter,
            }
        );
    }
    if ( $active_parameter_href->{frequency_filter} ) {

        my $infamily_directory = catdir( $active_parameter_href->{outdata_dir},
            $family_id, $outaligner_dir );
        my $outfamily_directory = $infamily_directory;

        ($xargs_file_counter) = analysis_frequency_filter_rio(
            {
                active_parameter_href   => $active_parameter_href,
                call_type               => $call_type,
                FILEHANDLE              => $FILEHANDLE,
                file_info_href          => $file_info_href,
                file_path               => $file_path,
                infamily_directory      => $infamily_directory,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                outfamily_directory     => $outfamily_directory,
                parameter_href          => $parameter_href,
                program_info_path       => $program_info_path,
                program_name            => q{frequency_filter},
                sample_info_href        => $sample_info_href,
                stderr_path        => $program_info_path . $DOT . q{stderr.txt},
                xargs_file_counter => $xargs_file_counter,
            }
        );
    }
    if ( $active_parameter_href->{varianteffectpredictor} ) {

        ($xargs_file_counter) = analysis_vep_rio(
            {
                active_parameter_href   => $active_parameter_href,
                call_type               => $call_type,
                FILEHANDLE              => $FILEHANDLE,
                file_info_href          => $file_info_href,
                file_path               => $file_path,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                parameter_href          => $parameter_href,
                program_info_path       => $program_info_path,
                program_name            => q{varianteffectpredictor},
                sample_info_href        => $sample_info_href,
                stderr_path        => $program_info_path . $DOT . q{stderr.txt},
                xargs_file_counter => $xargs_file_counter,
            }
        );
    }
    if ( $active_parameter_href->{vcfparser} ) {

        my $infamily_directory = catdir( $active_parameter_href->{outdata_dir},
            $family_id, $outaligner_dir );

        my $outfamily_directory = $infamily_directory;

        ($xargs_file_counter) = analysis_mip_vcfparser_rio(
            {
                active_parameter_href   => $active_parameter_href,
                call_type               => $call_type,
                FILEHANDLE              => $FILEHANDLE,
                file_info_href          => $file_info_href,
                file_path               => $file_path,
                infamily_directory      => $infamily_directory,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                outfamily_directory     => $outfamily_directory,
                parameter_href          => $parameter_href,
                program_info_path       => $program_info_path,
                program_name            => q{vcfparser},
                sample_info_href        => $sample_info_href,
                xargs_file_counter      => $xargs_file_counter,
            }
        );
    }
    if ( $active_parameter_href->{snpeff} ) {

        my $infamily_directory = catdir( $active_parameter_href->{outdata_dir},
            $family_id, $outaligner_dir );
        my $outfamily_directory = $infamily_directory;

        ($xargs_file_counter) = analysis_snpeff_rio(
            {
                active_parameter_href   => $active_parameter_href,
                call_type               => $call_type,
                FILEHANDLE              => $FILEHANDLE,
                file_info_href          => $file_info_href,
                file_path               => $file_path,
                job_id_href             => $job_id_href,
                infamily_directory      => $infamily_directory,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                outfamily_directory     => $outfamily_directory,
                parameter_href          => $parameter_href,
                program_info_path       => $program_info_path,
                program_name            => q{snpeff},
                sample_info_href        => $sample_info_href,
                xargs_file_counter      => $xargs_file_counter,
            }
        );
    }
    if ( $active_parameter_href->{rankvariant} ) {

        my $rank_program_name = q{rankvariant};

        if ( defined $parameter_href->{dynamic_parameter}{unaffected}
            && @{ $parameter_href->{dynamic_parameter}{unaffected} } eq
            @{ $active_parameter_href->{sample_ids} } )
        {

            $log->warn(
q{Only unaffected sample in pedigree - skipping genmod 'models', 'score' and 'compound'}
            );

            ($xargs_file_counter) = analysis_rankvariant_rio_unaffected(
                {
                    active_parameter_href   => $active_parameter_href,
                    call_type               => $call_type,
                    FILEHANDLE              => $FILEHANDLE,
                    file_info_href          => $file_info_href,
                    file_path               => $file_path,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    parameter_href          => $parameter_href,
                    program_info_path       => $program_info_path,
                    program_name            => $rank_program_name,
                    sample_info_href        => $sample_info_href,
                    xargs_file_counter      => $xargs_file_counter,
                }
            );
        }
        else {

            ($xargs_file_counter) = analysis_rankvariant_rio(
                {
                    active_parameter_href   => $active_parameter_href,
                    call_type               => $call_type,
                    FILEHANDLE              => $FILEHANDLE,
                    file_info_href          => $file_info_href,
                    file_path               => $file_path,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    job_id_href             => $job_id_href,
                    parameter_href          => $parameter_href,
                    program_info_path       => $program_info_path,
                    program_name            => $rank_program_name,
                    sample_info_href        => $sample_info_href,
                    xargs_file_counter      => $xargs_file_counter,
                }
            );
        }
    }
    if ( $active_parameter_href->{endvariantannotationblock} ) {

        ($xargs_file_counter) = analysis_endvariantannotationblock_rio(
            {
                active_parameter_href   => $active_parameter_href,
                FILEHANDLE              => $FILEHANDLE,
                call_type               => $call_type,
                file_info_href          => $file_info_href,
                file_path               => $file_path,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                parameter_href          => $parameter_href,
                program_info_path       => $program_info_path,
                program_name            => q{endvariantannotationblock},
                sample_info_href        => $sample_info_href,
                xargs_file_counter      => $xargs_file_counter,
            }
        );
    }
    return;
}

1;
