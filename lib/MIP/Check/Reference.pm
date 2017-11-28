package MIP::Check::Reference;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ fileparse };
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use List::MoreUtils qw { uniq };
use Readonly;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ check_bwa_prerequisites check_capture_file_prerequisites check_file_endings_to_build check_human_genome_file_endings check_human_genome_prerequisites check_if_processed_by_vt check_parameter_metafiles check_references_for_vt };
}

## Constants
Readonly my $DOT     => q{.};
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };
Readonly my $TAB     => qq{\t};

sub check_references_for_vt {

## Function : Check if vt has processed references
## Returns  : @to_process_references
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $log                   => Log object
##          : $parameter_href        => Parameter hash {REF}
##          : $vt_references_ref     => The references to check with vt {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $log;
    my $parameter_href;
    my $vt_references_ref;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        log            => { store => \$log },
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        vt_references_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$vt_references_ref,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Checked references
    my @checked_references;

    ## Store references to process later
    my @to_process_references;

    ## Avoid checking the same reference multiple times
    my %seen;

  PARAMETER_NAME:
    foreach my $parameter_name ( @{$vt_references_ref} ) {

      ASSOCIATED_PROGRAM:
        foreach my $associated_program (
            @{ $parameter_href->{$parameter_name}{associated_program} } )
        {

            ## Alias
            my $active_associated_program =
              $active_parameter_href->{$associated_program};

            next ASSOCIATED_PROGRAM if ( not $active_associated_program );

            ## If SCALAR data type
            if ( $parameter_href->{$parameter_name}{data_type} eq q{SCALAR} ) {

                my $annotation_file = $active_parameter_href->{$parameter_name};

                if ( not exists $seen{$annotation_file} ) {

                    ## Check if vt has processed references using regexp
                    @checked_references = check_if_processed_by_vt(
                        {
                            reference_file_path => $annotation_file,
                            log                 => $log,
                        }
                    );
                    push @to_process_references, @checked_references;
                }
                $seen{$annotation_file} = undef;
            }
            elsif ( $parameter_href->{$parameter_name}{data_type} eq q{ARRAY} )
            {
                ## ARRAY reference

              ANNOTION_FILE:
                foreach my $annotation_file (
                    @{ $active_parameter_href->{$parameter_name} } )
                {

                    if ( not exists $seen{$annotation_file} ) {

                        ## Check if vt has processed references using regexp
                        @checked_references = check_if_processed_by_vt(
                            {
                                reference_file_path => $annotation_file,
                                log                 => $log,
                            }
                        );
                    }
                    push @to_process_references, @checked_references;
                    $seen{$annotation_file} = undef;
                }
            }
            elsif ( $parameter_href->{$parameter_name}{data_type} eq q{HASH} ) {
                ## Hash reference

              ANNOTATION_FILE:
                for my $annotation_file (
                    keys %{ $active_parameter_href->{$parameter_name} } )
                {

                    if ( not exists $seen{$annotation_file} ) {

                        ## Check if vt has processed references using regexp
                        @checked_references = check_if_processed_by_vt(
                            {
                                reference_file_path => $annotation_file,
                                log                 => $log,
                            }
                        );
                    }
                    push @to_process_references, @checked_references;
                    $seen{$annotation_file} = undef;
                }
            }
        }
    }
    return uniq(@to_process_references);
}

sub check_if_processed_by_vt {

## Function : Check if vt has processed references using regexp
## Returns  : @process_references
## Arguments: $log                 => Log object
##          : $reference_file_path => The reference file path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $log;
    my $reference_file_path;

    my $tmpl = {
        log => {
            required => 1,
            defined  => 1,
            store    => \$log
        },
        reference_file_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$reference_file_path
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %vt_regexp;

    $vt_regexp{vt_decompose}{vcf_key} = q{OLD_MULTIALLELIC};
    $vt_regexp{vt_normalize}{vcf_key} = q{OLD_VARIANT};

    my @to_process_references;
    ## Downloaded and vt later (for downloadable references otherwise
    ## file existens error is thrown downstream)
    if ( not -e $reference_file_path ) {

        ## Do nothing since there is not ref file to check
        return;
    }

  VT_PARAMETER_NAME:
    foreach my $vt_parameter_name ( keys %vt_regexp ) {
        ## MIP flags

        ### Assemble perl regexp for detecting VT keys in vcf
        ## Execute perl
        my $regexp = q?perl -nae '?;

        ## Find vcf_key
        $regexp .=
          q?if($_=~/ID\=? . $vt_regexp{$vt_parameter_name}{vcf_key} . q?/) { ?;

        ## Write to stdout
        $regexp .= q?print $_} ?;

        ## If header is finished quit
        $regexp .= q?if($_=~/#CHROM/) {last}'?;

        ## Detect if vt program has processed reference
        my $ret = `bcftools view $reference_file_path | $regexp`;

        ## No trace of vt processing found
        if ( not $ret ) {

            ## Add reference for downstream processing
            push @to_process_references, $reference_file_path;
            $log->warn( $TAB
                  . q{Cannot detect that }
                  . $vt_parameter_name
                  . q{ has processed reference: }
                  . $reference_file_path
                  . $NEWLINE );
        }
        else {
            ## Found vt processing trace

            $log->info( $TAB
                  . q{Reference check: }
                  . $reference_file_path
                  . q{ vt: }
                  . $vt_parameter_name
                  . q{ - PASS}
                  . $NEWLINE );
        }
    }
    return uniq(@to_process_references);
}

sub check_human_genome_prerequisites {

## Function : Checks if the human genome prerequisites needs to be built and builds them if required
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $log                     => Log object
##          : $program_name            => Program name
##          : $parameter_href          => Parameter hash {REF}
##          : $sample_info_href        => Info on samples and family hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $log;
    my $parameter_href;
    my $program_name;
    my $sample_info_href;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
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
        log => {
            required => 1,
            defined  => 1,
            store    => \$log
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Recipes::Build::Human_genome_prerequisites
      qw{ build_human_genome_prerequisites };

    ## Files assocaiated with human genome reference
    if (   $parameter_href->{human_genome_reference}{build_file} == 1
        || $file_info_href->{human_genome_compressed} eq q{compressed} )
    {

        ## Creates the humanGenomePreRequisites using active_parameters{human_genome_reference} as reference.
        build_human_genome_prerequisites(
            {
                parameter_href          => $parameter_href,
                active_parameter_href   => $active_parameter_href,
                sample_info_href        => $sample_info_href,
                file_info_href          => $file_info_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                family_id               => $active_parameter_href->{family_id},
                outaligner_dir => $active_parameter_href->{outaligner_dir},
                program_name   => $program_name,
                log            => $log,
            }
        );
    }
    return 1;
}

sub check_human_genome_file_endings {

## Function : Check the existance of associated human genome files.
## Returns  :
## Arguments: $active_parameter_href              => Holds all set parameter for analysis {REF}
##          : $file_info_href                     => File info hash {REF}
##          : $human_genome_reference_name_prefix => The associated human genome file without file ending {REF}
##          : $parameter_href                     => Parameter hash {REF}
##          : $parameter_name                     => The parameter under evaluation
##          : $reference_dir                      => MIP reference directory {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $parameter_href;
    my $parameter_name;

    ## Default(s)
    my $human_genome_reference_name_prefix;
    my $reference_dir;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        human_genome_reference_name_prefix => {
            default =>
              $arg_href->{file_info_href}{human_genome_reference_name_prefix},
            strict_type => 1,
            store       => \$human_genome_reference_name_prefix
        },
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        parameter_name => { strict_type => 1, store => \$parameter_name },
        reference_dir  => {
            default     => $arg_href->{active_parameter_href}{reference_dir},
            strict_type => 1,
            store       => \$reference_dir,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Path qw{ check_filesystem_objects_existance };
    use MIP::Get::File qw{ get_seq_dict_contigs };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    ## Count the number of files that exists
    my $existence_check_counter = 0;

    ## Unpack
    my @human_reference_file_endings =
      @{ $file_info_href->{human_genome_reference_file_endings} };

  FILE_ENDING:
    foreach my $file_ending (@human_reference_file_endings) {

        my $path = $active_parameter_href->{human_genome_reference};

        if ( $file_ending eq q{.dict} ) {

            ## Removes ".file_ending" in filename.FILENDING(.gz)
            my ( $file, $dir_path ) =
              fileparse( $path, qr/ [.]fasta | [.]fasta[.]gz /sxm );
            $path = catfile( $dir_path, $file );
        }

        #Add current ending
        $path = $path . $file_ending;

        my ($does_exist) = check_filesystem_objects_existance(
            {
                object_name    => $path,
                parameter_name => $parameter_name,
                object_type    => q{file},
            }
        );

        ## Sum up the number of file that exists
        $existence_check_counter = $existence_check_counter + $does_exist;
    }
    ## Files need to be built
    if ( $existence_check_counter != scalar @human_reference_file_endings ) {

        $parameter_href->{$parameter_name}{build_file} = 1;
    }
    else {

        # All file exist in this check
        $parameter_href->{$parameter_name}{build_file} = 0;

        ## Get sequence contigs from human reference ".dict" file since it exists
        my $dict_file_path = catfile( $reference_dir,
            $human_genome_reference_name_prefix . $DOT . q{dict} );
        @{ $file_info_href->{contigs} } = get_seq_dict_contigs(
            {
                dict_file_path => $dict_file_path,
                log            => $log,
            }
        );
    }
    return;
}

sub check_capture_file_prerequisites {

## Function : Check if capture file prerequisites needs to be built and initate process to build them if required
## Returns  :
## Arguments: $active_parameter_href       => Active parameters for this analysis hash {REF}
##          : $FILEHANDLE                  => Filehandle to write to
##          : $infile_lane_prefix_href     => Infile(s) without the ".ending" {REF}
##          : $infile_list_suffix          => Infile list suffix
##          : $job_id_href                 => Job id hash {REF}
##          : $log                         => Log object
##          : $padded_infile_list_suffix   => Padded infile list suffix
##          : $padded_interval_list_suffix => Padded interval list suffix
##          : $parameter_href              => Parameter hash {REF}
##          : $program_name                => Program name
##          : $sample_info_href            => Info on samples and family hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $FILEHANDLE;
    my $infile_lane_prefix_href;
    my $infile_list_suffix;
    my $job_id_href;
    my $log;
    my $padded_infile_list_suffix;
    my $padded_interval_list_suffix;
    my $parameter_href;
    my $program_name;
    my $sample_info_href;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        FILEHANDLE              => { store => \$FILEHANDLE, },
        infile_lane_prefix_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$infile_lane_prefix_href,
        },
        infile_list_suffix => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_list_suffix,
        },
        job_id_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$job_id_href,
        },
        log => {
            required => 1,
            defined  => 1,
            store    => \$log,
        },
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        padded_infile_list_suffix => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$padded_infile_list_suffix,
        },
        padded_interval_list_suffix => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$padded_interval_list_suffix,
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Recipes::Build::Capture_file_prerequisites
      qw{ build_capture_file_prerequisites };

    if ( $parameter_href->{exome_target_bed}{build_file} == 1 ) {

        build_capture_file_prerequisites(
            {
                parameter_href              => $parameter_href,
                active_parameter_href       => $active_parameter_href,
                sample_info_href            => $sample_info_href,
                infile_lane_prefix_href     => $infile_lane_prefix_href,
                job_id_href                 => $job_id_href,
                infile_list_suffix          => $infile_list_suffix,
                padded_infile_list_suffix   => $padded_infile_list_suffix,
                padded_interval_list_suffix => $padded_interval_list_suffix,
                program_name                => $program_name,
                FILEHANDLE                  => $FILEHANDLE,
                log                         => $log,
            }
        );

        ## Only build once for all modules and files
        $parameter_href->{exome_target_bed}{build_file} = 0;
    }
    return;
}

sub check_bwa_prerequisites {

## Function : Checks if the Bwa mem prerequisites needs to be built and builds them if required
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $parameter_build_name          => File info build key parameter name
##          : $parameter_href          => Parameter hash {REF}
##          : $program_name            => Program name
##          : $sample_info_href        => Info on samples and family hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $parameter_build_name;
    my $parameter_href;
    my $program_name;
    my $sample_info_href;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
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
        parameter_build_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$parameter_build_name,
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Recipes::Build::Bwa_prerequisites qw{ build_bwa_prerequisites };

    if ( $parameter_href->{bwa_build_reference}{build_file} == 1 ) {

        build_bwa_prerequisites(
            {
                parameter_href          => $parameter_href,
                active_parameter_href   => $active_parameter_href,
                sample_info_href        => $sample_info_href,
                file_info_href          => $file_info_href,
                infile_lane_prefix_href => $infile_lane_prefix_href,
                job_id_href             => $job_id_href,
                bwa_build_reference_file_endings_ref =>
                  \@{ $file_info_href->{$parameter_build_name} },
                program_name => $program_name,
            }
        );
    }
    return;
}

sub check_file_endings_to_build {

## Function : Checks files to be built by combining filename stub with fileendings.
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $file_endings_ref      => Reference to the file_endings to be added to the filename stub {REF}
##          : $file_name             => File name
##          : $parameter_href        => Parameter hash {REF}
##          : $parameter_name        => MIP parameter name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_endings_ref;
    my $file_name;
    my $parameter_href;
    my $parameter_name;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        file_endings_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$file_endings_ref
        },
        file_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$file_name
        },
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
        parameter_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$parameter_name,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Check::Path qw{ check_filesystem_objects_existance };

    ## Count the number of files that exists
    my $existence_check_counter = 0;

  FILE_ENDING:
    foreach my $file_ending ( @{$file_endings_ref} ) {

        my ($exist) = check_filesystem_objects_existance(
            {
                object_name    => catfile( $file_name . $file_ending ),
                parameter_name => $parameter_name,
                object_type    => q{file},
            }
        );

        ## Sum up the number of file that exists
        $existence_check_counter = $existence_check_counter + $exist;
    }

    ## Files need to be built
    if ( $existence_check_counter != scalar @{$file_endings_ref} ) {

        $parameter_href->{$parameter_name}{build_file} = 1;
    }
    else {

        # All file exist in this check
        $parameter_href->{$parameter_name}{build_file} = 0;
    }
    return;
}

sub check_parameter_metafiles {

## Function : Checks parameter metafile exists
## Returns  :
## Arguments: $active_parameter_href => Holds all set parameter for analysis
##          : $file_info_href        => File info hash {REF}
##          : $parameter_href        => Holds all parameters

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $parameter_href;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href,
        },
        parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$parameter_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  PARAMETER:
    foreach my $parameter_name ( %{$file_info_href} ) {

        ## Active parameter
        my $parameter = $active_parameter_href->{$parameter_name};

        next PARAMETER if ( not $parameter );

      ASSOCIATED_PROGRAM:
        foreach my $associated_program (
            @{ $parameter_href->{$parameter_name}{associated_program} } )
        {

            ## Active associated program
            my $active_associated_program =
              $active_parameter_href->{$associated_program};

            ## Only check active parmeters
            next ASSOCIATED_PROGRAM
              if ( not defined $active_associated_program );

            if ( ref $parameter eq q{HASH} ) {

              PATH:
                for my $path ( keys %{$parameter} ) {

                    ## Checks files to be built by combining filename stub with fileendings
                    check_file_endings_to_build(
                        {
                            parameter_href        => $parameter_href,
                            active_parameter_href => $active_parameter_href,
                            file_endings_ref =>
                              \@{ $file_info_href->{$parameter_name} },
                            parameter_name => $parameter_name,
                            file_name      => $path,
                        }
                    );
                }
            }
            else {

                ## Checks files to be built by combining filename stub with fileendings
                check_file_endings_to_build(
                    {
                        parameter_href        => $parameter_href,
                        active_parameter_href => $active_parameter_href,
                        file_endings_ref =>
                          \@{ $file_info_href->{$parameter_name} },
                        parameter_name => $parameter_name,
                        file_name =>
                          $active_parameter_href->{human_genome_reference},
                    }
                );
            }
        }
    }
    return;
}

1;
