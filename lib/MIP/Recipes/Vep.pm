package MIP::Recipes::Vep;

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
use File::Basename qw{ fileparse };
use File::Spec::Functions qw{ catdir catfile };

## CPANM
use Readonly;

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw(analysis_vep);

}

##Constants
Readonly my $NEWLINE    => qq{\n};
Readonly my $UNDERSCORE => q{_};

sub analysis_vep {

##analysis_vep

##Function : Varianteffectpredictor performs effect predictions and annotation of variants.
##Returns  : "|$xargs_file_counter"
##Arguments: $parameter_href, $active_parameter_href, $sample_info_href, $file_info_href, $infile_lane_prefix_href, $job_id_href, $program_name, $program_info_path, $file_path, $stderr_path, $FILEHANDLE, family_id_ref, $temp_directory_ref, $outaligner_dir_ref, $call_type, $xargs_file_counter
##         : $parameter_href          => Parameter hash {REF}
##         : $active_parameter_href   => Active parameters for this analysis hash {REF}
##         : $sample_info_href        => Info on samples and family hash {REF}
##         : $file_info_href          => The file_info hash {REF}
##         : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##         : $job_id_href             => Job id hash {REF}
##         : $program_name            => Program name
##         : $program_info_path       => The program info path
##         : $file_path               => File path
##         : $stderr_path             => Stderr path of the block script
##         : $FILEHANDLE              => Filehandle to write to
##         : $family_id_ref           => Family id {REF}
##         : $temp_directory_ref      => Temporary directory {REF}
##         : $outaligner_dir_ref      => Outaligner_dir used in the analysis {REF}
##         : $call_type               => The variant call type
##         : $xargs_file_counter      => The xargs file counter

    my ($arg_href) = @_;

    ## Default(s)
    my $family_id_ref;
    my $temp_directory_ref;
    my $outaligner_dir_ref;
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
    my $program_info_path;
    my $file_path;
    my $stderr_path;
    my $FILEHANDLE;

    my $tmpl = {
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href
        },
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href
        },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href
        },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href
        },
        program_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$program_name
        },
        program_info_path => { strict_type => 1, store => \$program_info_path },
        file_path         => { strict_type => 1, store => \$file_path },
        stderr_path       => { strict_type => 1, store => \$stderr_path },
        FILEHANDLE    => { store => \$FILEHANDLE },
        family_id_ref => {
            default     => \$arg_href->{active_parameter_href}{family_id},
            strict_type => 1,
            store       => \$family_id_ref
        },
        temp_directory_ref => {
            default     => \$arg_href->{active_parameter_href}{temp_directory},
            strict_type => 1,
            store       => \$temp_directory_ref
        },
        outaligner_dir_ref => {
            default     => \$arg_href->{active_parameter_href}{outaligner_dir},
            strict_type => 1,
            store       => \$outaligner_dir_ref
        },
        call_type =>
          { default => "BOTH", strict_type => 1, store => \$call_type },
        xargs_file_counter => {
            default     => 0,
            allow       => qr/ ^\d+$ /xsm,
            strict_type => 1,
            store       => \$xargs_file_counter
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Cluster qw(get_core_number);
    use MIP::IO::Files qw(migrate_file xargs_migrate_contig_files);
    use MIP::Set::File qw{set_file_suffix};
    use MIP::Get::File qw{get_file_suffix};
    use MIP::Recipes::Xargs qw{ xargs_command };
    use Program::Variantcalling::Vep qw(variant_effect_predictor);
    use MIP::QC::Record qw(add_program_outfile_to_sample_info);
    use MIP::Processmanagement::Slurm_processes
      qw(slurm_submit_job_sample_id_dependency_add_to_family);

    my $reduce_io_ref = \$active_parameter_href->{reduce_io};
    my $xargs_file_path_prefix;
    my $job_id_chain = $parameter_href->{ "p" . $program_name }{chain};

    ## Filehandles
    my $XARGSFILEHANDLE = IO::Handle->new();    #Create anonymous filehandle

    unless ( defined($FILEHANDLE) ) {           #Run as individual sbatch script

        $FILEHANDLE = IO::Handle->new();        #Create anonymous filehandle
    }

    ## Get core number depending on user supplied input exists or not and max number of cores
    my $core_number = get_core_number(
        {
            module_core_number => $active_parameter_href->{module_core_number}
              { "p" . $program_name },
            modifier_core_number => scalar( @{ $file_info_href->{contigs} } ),
            max_cores_per_node => $active_parameter_href->{max_cores_per_node},
        }
    );

    my $fork_number = 4;    #Varianteffectpredictor forks
    $core_number =
      floor( $core_number / $fork_number );    #Adjust for the number of forks

    if ( !$$reduce_io_ref ) {                  #Run as individual sbatch script

        use MIP::Script::Setup_script qw(setup_script);

        ## Creates program directories (info & programData & programScript), program script filenames and writes sbatch header
        ( $file_path, $program_info_path ) = setup_script(
            {
                active_parameter_href => $active_parameter_href,
                job_id_href           => $job_id_href,
                FILEHANDLE            => $FILEHANDLE,
                directory_id          => $$family_id_ref,
                program_name          => $program_name,
                program_directory     => catfile( lc($$outaligner_dir_ref) ),
                call_type             => $call_type,
                core_number           => $core_number,
                process_time =>
                  $active_parameter_href->{module_time}{ "p" . $program_name },
                temp_directory => $$temp_directory_ref
            }
        );
        $stderr_path = $program_info_path . ".stderr.txt";
    }

    # Split to enable submission to &sample_info_qc later
    my ( $volume, $directory, $stderr_file ) =
      File::Spec->splitpath($stderr_path);

    ## Assign directories
    my $infamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref );
    my $outfamily_directory = catdir( $active_parameter_href->{outdata_dir},
        $$family_id_ref, $$outaligner_dir_ref );
    $parameter_href->{ "p" . $program_name }{indirectory} =
      $outfamily_directory;    #Used downstream

    ## Assign file_tags
    my $infile_tag = $file_info_href->{$$family_id_ref}{pvt}{file_tag};
    my $outfile_tag =
      $file_info_href->{$$family_id_ref}{ "p" . $program_name }{file_tag};
    my $infile_prefix       = $$family_id_ref . $infile_tag . $call_type;
    my $file_path_prefix    = catfile( $$temp_directory_ref, $infile_prefix );
    my $outfile_prefix      = $$family_id_ref . $outfile_tag . $call_type;
    my $outfile_path_prefix = catfile( $$temp_directory_ref, $outfile_prefix );

    ### Assign suffix
    ## Return the current infile vcf compression suffix for this jobid chain
    my $infile_suffix = get_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_file_suffix",
            jobid_chain    => $job_id_chain,
        }
    );
    my $outfile_suffix = set_file_suffix(
        {
            parameter_href => $parameter_href,
            suffix_key     => "variant_file_suffix",
            job_id_chain   => $job_id_chain,
            file_suffix =>
              $parameter_href->{ "p" . $program_name }{outfile_suffix},
        }
    );

    if ( !$$reduce_io_ref ) {    #Run as individual sbatch script

        ## Copy file(s) to temporary directory
        say $FILEHANDLE "## Copy file(s) to temporary directory";
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
                indirectory        => $infamily_directory,
                temp_directory     => $$temp_directory_ref,
            }
        );
    }

    ## varianteffectpredictor
    say $FILEHANDLE "## Varianteffectpredictor";

    my $assembly_version = $file_info_href->{human_genome_reference_source}
      . $file_info_href->{human_genome_reference_version};

    ## Alias genome source and version to be compatible with VEP
    alias_assembly_version( { assembly_version_ref => \$assembly_version } );

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

    foreach my $contig ( @{ $file_info_href->{contigs_size_ordered} } ) {

        ## Get parameters
        ## VEP plugins
        my @plugins;
        foreach my $plugin ( @{ $active_parameter_href->{vep_plugins} } ) {

            if ( $plugin eq "LoF" ) {

                push(
                    @plugins,
                    $plugin
                      . ",human_ancestor_fa:"
                      . catfile(
                        $active_parameter_href->{vep_directory_cache},
                        "human_ancestor.fa,filter_position:0.05"
                      )
                );
            }
            elsif ( $plugin eq "UpDownDistance" )
            {    #Special case for mitochondrial contig annotation

                if ( $contig =~ /MT|M/ ) {

                    push( @plugins, "UpDownDistance,10,10" );
                }
            }
            else {

                push( @plugins, $plugin );
            }
        }

        ## VEPFeatures
        my @vep_features_ref;
        foreach my $vep_feature ( @{ $active_parameter_href->{vep_features} } )
        {

            push( @vep_features_ref, $vep_feature )
              ;    #Add VEP features to the output.

            if ( ( $contig =~ /MT|M/ ) && ( $vep_feature eq "refseq" ) )
            {      #Special case for mitochondrial contig annotation

                push( @vep_features_ref, "all_refseq" );
            }
        }

        variant_effect_predictor(
            {
                regions_ref      => [$contig],
                plugins_ref      => \@plugins,
                vep_features_ref => \@vep_features_ref,
                script_path      => catfile(
                    $active_parameter_href->{vep_directory_path},
                    "variant_effect_predictor.pl"
                ),
                assembly => $assembly_version,
                cache_directory =>
                  $active_parameter_href->{vep_directory_cache},
                reference_path =>
                  $active_parameter_href->{human_genome_reference},
                infile_format  => substr( $outfile_suffix, 1 ),
                outfile_format => substr( $outfile_suffix, 1 ),
                fork           => $fork_number,
                buffer_size    => 20000,
                infile_path    => $file_path_prefix . "_"
                  . $contig
                  . $infile_suffix,
                outfile_path => $outfile_path_prefix . "_"
                  . $contig
                  . $infile_suffix,
                stderrfile_path => $xargs_file_path_prefix . "."
                  . $contig
                  . ".stderr.txt",
                stdoutfile_path => $xargs_file_path_prefix . "."
                  . $contig
                  . ".stdout.txt",
                FILEHANDLE => $XARGSFILEHANDLE,
            }
        );
        print $XARGSFILEHANDLE "\n";
    }

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        ## Collect QC metadata info for later use
        my $qc_vep_summary_outfile =
            $outfile_prefix . q{_}
          . $file_info_href->{contigs_size_ordered}[0]
          . $infile_suffix
          . q{_summary.html};
        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => $program_name . q{summary},
                outdirectory     => $outfamily_directory,
                outfile          => $qc_vep_summary_outfile,
            }
        );
        ## Collect QC metadata info for later use
        my $qc_vep_outfile =
            $outfile_prefix . q{_}
          . $file_info_href->{contigs_size_ordered}[0]
          . $infile_suffix;
        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => $program_name,
                outdirectory     => $outfamily_directory,
                outfile          => $qc_vep_outfile,
            }
        );
    }

    ## QC Data File(s)
    migrate_file(
        {
            infile_path => $outfile_path_prefix . q{_*}
              . $infile_suffix . q{_s*},
            outfile_path => $outfamily_directory,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say $FILEHANDLE q{wait}, "\n";

    close($XARGSFILEHANDLE);

    if ( !$$reduce_io_ref ) {    #Run as individual sbatch script

        ## Copies file from temporary directory.
        say $FILEHANDLE "## Copy file from temporary directory";
        migrate_file(
            {
                infile_path => $outfile_path_prefix . q{_*}
                  . $infile_suffix . q{*},
                outfile_path => $outfamily_directory,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
        say $FILEHANDLE q{wait}, "\n";

        close($FILEHANDLE);
    }
    else {    #Move file for downstream collection of VEP version

        ## Copies file from temporary directory.
        say $FILEHANDLE "## Copy file from temporary directory";
        migrate_file(
            {
                infile_path => $outfile_path_prefix . q{_}
                  . $file_info_href->{contigs_size_ordered}[0]
                  . $infile_suffix,
                outfile_path => $outfamily_directory,
                FILEHANDLE   => $FILEHANDLE,
            }
        );
        say $FILEHANDLE q{wait}, "\n";
    }

    if ( $active_parameter_href->{ "p" . $program_name } == 1 ) {

        if ( !$$reduce_io_ref ) {    #Run as individual sbatch script

            slurm_submit_job_sample_id_dependency_add_to_family(
                {
                    job_id_href             => $job_id_href,
                    infile_lane_prefix_href => $infile_lane_prefix_href,
                    sample_ids_ref =>
                      \@{ $active_parameter_href->{sample_ids} },
                    family_id        => $$family_id_ref,
                    path             => $job_id_chain,
                    log              => $log,
                    sbatch_file_name => $file_path,
                }
            );
        }
        if ($$reduce_io_ref)
        { #Redirect qccollect search to Block File, since VEP will write stderr there

            $program_name = "variantannotationblock";
        }

        ## Collect QC metadata info for later use
        add_program_outfile_to_sample_info(
            {
                sample_info_href => $sample_info_href,
                program_name     => $program_name,
                outdirectory     => $directory,
                outfile          => $stderr_file,
            }
        );
    }
    if ($$reduce_io_ref) {

        return
          $xargs_file_counter
          ; #Track the number of created xargs scripts per module for Block algorithm
    }
}

1;
