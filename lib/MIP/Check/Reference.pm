package MIP::Check::Reference;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };

## CPANM
use Readonly;
use List::MoreUtils qw { uniq };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ check_references_for_vt check_if_processed_by_vt check_human_genome_prerequisites check_capture_file_prerequisites check_bwa_prerequisites };
}

## Constants
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE   => q{ };
Readonly my $TAB     => qq{\t};

sub check_references_for_vt {

## Function : Check if vt has processed references
## Returns  : @to_process_references
## Arguments: $parameter_href, $active_parameter_href, $vt_references_ref
##          : $parameter_href        => Parameter hash {REF}
##          : $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $vt_references_ref     => The references to check with vt {REF}
##          : $log                   => Log object

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $vt_references_ref;
    my $log;

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
        vt_references_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$vt_references_ref,
        },
        log => { store => \$log },
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

##check_if_processed_by_vt

##Function : Check if vt has processed references using regexp
##Returns  : @process_references
##Arguments: $reference_file_path, $log
##         : $reference_file_path => The reference file path
##         : $log                 => Log object
    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $reference_file_path;
    my $log;

    my $tmpl = {
        reference_file_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$reference_file_path
        },
        log => {
            required => 1,
            defined  => 1,
            store    => \$log
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
        my $ret = `less $reference_file_path | $regexp`;

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
## Arguments: $parameter_href          => Parameter hash {REF}
##          : $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $program_name            => Program name
##          : $log                     => Log object

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $file_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $program_name;
    my $log;

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
        log => {
            required => 1,
            defined  => 1,
            store    => \$log
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

sub check_capture_file_prerequisites {

## Function : Check if capture file prerequisites needs to be built and initate process to build them if required
## Returns  :
## Arguments: $parameter_href              => Parameter hash {REF}
##          : $active_parameter_href       => Active parameters for this analysis hash {REF}
##          : $sample_info_href            => Info on samples and family hash {REF}
##          : $infile_lane_prefix_href     => Infile(s) without the ".ending" {REF}
##          : $job_id_href                 => Job id hash {REF}
##          : $infile_list_suffix          => Infile list suffix
##          : $padded_infile_list_suffix   => Padded infile list suffix
##          : $padded_interval_list_suffix => Padded interval list suffix
##          : $program_name                => Program name
##          : $FILEHANDLE                  => Filehandle to write to
##          : $log                         => Log object

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;
    my $active_parameter_href;
    my $sample_info_href;
    my $infile_lane_prefix_href;
    my $job_id_href;
    my $infile_list_suffix;
    my $padded_infile_list_suffix;
    my $padded_interval_list_suffix;
    my $program_name;
    my $FILEHANDLE;
    my $log;

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
        infile_list_suffix => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_list_suffix,
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
        FILEHANDLE => { store => \$FILEHANDLE, },
        log        => {
            required => 1,
            defined  => 1,
            store    => \$log,
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
## Arguments: $parameter_href          => Parameter hash {REF}
##          : $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $sample_info_href        => Info on samples and family hash {REF}
##          : $file_info_href          => File info hash {REF}
##          : $infile_lane_prefix_href => Infile(s) without the ".ending" {REF}
##          : $job_id_href             => Job id hash {REF}
##          : $program_name            => Program name

    my ($arg_href) = @_;

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
                  \@{ $file_info_href->{bwa_build_reference_file_endings} },
                program_name => $program_name,
            }
        );
    }
    return;
}

1;
