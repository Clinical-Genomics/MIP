package MIP::Check::Reference;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use List::MoreUtils qw { uniq };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $LOG_NAME $NEWLINE $SPACE $TAB };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      check_if_processed_by_bcftools
      check_references_for_bcftools
    };
}

sub check_if_processed_by_bcftools {

## Function : Check if bcftools has processed references using regexp
## Returns  : @process_references
## Arguments: $reference_file_path  => The reference file path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $reference_file_path;

    my $tmpl = {
        reference_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$reference_file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Format::Vcf qw{ get_bcftools_norm_command_from_vcf_header };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my @to_process_references;
    ## Downloaded and check bcftools header later (for downloadable references otherwise
    ## file existens error is thrown downstream)
    return if ( not -e $reference_file_path );

    my $header_id_line = get_bcftools_norm_command_from_vcf_header(
        {
            vcf_file_path => $reference_file_path,
        }
    );

    ## No trace of bcftools processing found
    if ( not $header_id_line ) {

        ## Add reference for downstream processing
        push @to_process_references, $reference_file_path;
        $log->warn( $TAB
              . q{Cannot detect that bcftools norm}
              . q{ has processed reference: }
              . $reference_file_path
              . $NEWLINE );
    }
    else {

        ## Found bcftools processing trace
        $log->info( $TAB
              . q{Reference check: }
              . $reference_file_path
              . q{ bcftools norm: }
              . q{ - PASS}
              . $NEWLINE );
    }

    return uniq(@to_process_references);
}

sub check_references_for_bcftools {

## Function : Check if bcftools has processed references
## Returns  : @to_process_references
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $bcftools_references_ref => The references to check with bcftools {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $parameter_href;
    my $bcftools_references_ref;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        bcftools_references_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$bcftools_references_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Parameter qw{ get_parameter_attribute };

    ## Checked references
    my @checked_references;

    ## Store references to process later
    my @to_process_references;

    ## Avoid checking the same reference multiple times
    my %seen;

    ## TOML parameters
    my %toml = (
        sv_vcfanno_config => 1,
        vcfanno_config    => 1,
    );

  PARAMETER_NAME:
    foreach my $parameter_name ( @{$bcftools_references_ref} ) {

        my @attributes = get_parameter_attribute(
            {
                attribute      => q{associated_recipe},
                parameter_href => $parameter_href,
                parameter_name => $parameter_name,
            }
        );

      ASSOCIATED_RECIPE:
        foreach my $associated_recipe (@attributes) {

            ## Alias
            my $active_associated_recipe = $active_parameter_href->{$associated_recipe};

            next ASSOCIATED_RECIPE if ( not $active_associated_recipe );

            ## If SCALAR data type
            if ( $parameter_href->{$parameter_name}{data_type} eq q{SCALAR} ) {

                my $annotation_file = $active_parameter_href->{$parameter_name};

                ## Special case for toml configs (annotation file path recorded inside file parameter)
                if ( defined $toml{$parameter_name} ) {

                    _parse_vcfanno_toml_path(
                        {
                            seen_href                 => \%seen,
                            toml_file_path            => $annotation_file,
                            to_process_references_ref => \@to_process_references,
                        }
                    );
                }
                if ( not exists $seen{$annotation_file} ) {

                    ## Check if bcftools has processed references using regexp
                    @checked_references = check_if_processed_by_bcftools(
                        {
                            reference_file_path => $annotation_file,
                        }
                    );
                    push @to_process_references, @checked_references;
                }
                $seen{$annotation_file} = undef;
            }
            elsif ( $parameter_href->{$parameter_name}{data_type} eq q{ARRAY} ) {
                ## ARRAY reference

              ANNOTION_FILE:
                foreach my $annotation_file ( @{ $active_parameter_href->{$parameter_name} } ) {

                    if ( not exists $seen{$annotation_file} ) {

                        ## Check if bcftools has processed references using regexp
                        @checked_references = check_if_processed_by_bcftools(
                            {
                                reference_file_path => $annotation_file,
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
                for my $annotation_file ( keys %{ $active_parameter_href->{$parameter_name} } ) {

                    if ( not exists $seen{$annotation_file} ) {

                        ## Check if bcftools has processed references using regexp
                        @checked_references = check_if_processed_by_bcftools(
                            {
                                reference_file_path => $annotation_file,
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

sub _parse_vcfanno_toml_path {

## Function : Parse TOML config for path to check with bcftools
## Returns  :
## Arguments: $seen_href                 => Avoid checking the same reference multiple times
##          : $toml_file_path            => Toml config file path
##          : $to_process_references_ref => Store references to process later

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $seen_href;
    my $toml_file_path;
    my $to_process_references_ref;

    my $tmpl = {
        seen_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$seen_href,
            strict_type => 1,
        },
        to_process_references_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$to_process_references_ref,
            strict_type => 1,
        },
        toml_file_path => {
            default     => 1,
            store       => \$toml_file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Io::Read qw{ read_from_file };

    my %vcfanno_config = read_from_file(
        {
            format => q{toml},
            path   => $toml_file_path,
        }
    );

    ## Add config parameter to avoid bcftools check of toml config path
    $seen_href->{$toml_file_path} = undef;

  ANNOTATION:
    foreach my $annotation_href ( @{ $vcfanno_config{annotation} } ) {

        ## Annotation file path to check
        my $annotation_file_path = $annotation_href->{file};

        ## Only check vcf files which should have the fields annotation
        next ANNOTATION if ( not exists $annotation_href->{fields} );

        ## Skip files not set for normalization in the toml
        next ANNOTATION if ( $annotation_href->{skip_split_and_normalize} );

        if ( not exists $seen_href->{$annotation_file_path} ) {

            ## Check if bcftools has processed references using regexp
            my @checked_references = check_if_processed_by_bcftools(
                {
                    reference_file_path => $annotation_file_path,
                }
            );
            push @{$to_process_references_ref}, @checked_references;
        }
        $seen_href->{$annotation_file_path} = undef;
    }
    return;
}

1;
