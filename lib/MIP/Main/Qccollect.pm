package MIP::Main::Qccollect;

#### Collects MPS QC from MIP. Loads information on files to examine and values
#### to extract from in YAML format and outputs exracted metrics in YAML format.

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use File::Basename qw{ fileparse };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ open close :all };
use Modern::Perl qw{ 2018 };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $LOG_NAME };
use MIP::Io::Read qw{ read_from_file };
use MIP::Io::Write qw{ write_to_file };
use MIP::Qccollect qw{ evaluate_analysis };
use MIP::Qc_data qw{ set_qc_data_recipe_info };
use MIP::Store qw{ store_metrics };

BEGIN {

    use base qw{ Exporter };
    require Exporter;

    # Set the version for version checking
    our $VERSION = q{2.1.10};

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ mip_qccollect parse_limit_qc_output };
}

sub mip_qccollect {

## Function : Execute mip qccollect of data metrics
## Returns  :
## Arguments: $eval_metric_file      => File with evaluation metrics
##          : $evaluate_plink_gender => Evaluate plink gender
##          : $limit_qc_output       => Only print a limited number of qc metrics
##          : $outfile               => Data metric output file
##          : $regexp_file           => Regular expression file
##          : $sample_info_file      => Sample info file
##          : $skip_evaluation       => Skip evaluation step
##          : $store_metrics_outfile => Store metrics outfile

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $eval_metric_file;
    my $evaluate_plink_gender;
    my $limit_qc_output;
    my $outfile;
    my $regexp_file;
    my $sample_info_file;
    my $skip_evaluation;
    my $store_metrics_outfile;

    my $tmpl = {
        eval_metric_file => {
            store       => \$eval_metric_file,
            strict_type => 1,
        },
        evaluate_plink_gender => {
            store       => \$evaluate_plink_gender,
            strict_type => 1,
        },
        outfile => {
            defined     => 1,
            required    => 1,
            store       => \$outfile,
            strict_type => 1,
        },
        limit_qc_output => {
            store       => \$limit_qc_output,
            strict_type => 1,
        },
        regexp_file => {
            defined     => 1,
            required    => 1,
            store       => \$regexp_file,
            strict_type => 1,
        },
        sample_info_file => {
            defined     => 1,
            required    => 1,
            store       => \$sample_info_file,
            strict_type => 1,
        },
        skip_evaluation => {
            store       => \$skip_evaluation,
            strict_type => 1,
        },
        store_metrics_outfile => {
            store       => \$store_metrics_outfile,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Save final output data
    my %qc_data;

    ## Save header(s) in each outfile
    my %qc_header;

    ## Save data in each outfile
    my %qc_recipe_data;

    ## Loads a YAML file into an arbitrary hash and returns it
    $log->info( q{Loading: } . $sample_info_file );
    my %sample_info = read_from_file(
        {
            format => q{yaml},
            path   => $sample_info_file,
        }
    );
    $log->info( q{Loaded: } . $sample_info_file );

    ## Loads a reg exp file into an arbitrary hash
    $log->info( q{Loading: } . $regexp_file );
    my %regexp = read_from_file(
        {
            format => q{yaml},
            path   => $regexp_file,
        }
    );
    $log->info( q{Loaded: } . $regexp_file );

    ## Set qccollect version to qc_data hash
    set_qc_data_recipe_info(
        {
            key          => q{version},
            qc_data_href => \%qc_data,
            recipe_name  => q{qccollect},
            value        => $MIP::Main::Qccollect::VERSION,
        }
    );

    ## Set supplied regexp file to qc_data hash
    set_qc_data_recipe_info(
        {
            key          => q{regexp_file},
            qc_data_href => \%qc_data,
            recipe_name  => q{qccollect},
            value        => $regexp_file,
        }
    );

    ## Extracts all Qc data on sample_id level using information in %sample_info and %regexp
    sample_qc(
        {
            qc_data_href        => \%qc_data,
            qc_header_href      => \%qc_header,
            qc_recipe_data_href => \%qc_recipe_data,
            regexp_href         => \%regexp,
            sample_info_href    => \%sample_info,
        }
    );

    ## Extracts all Qc data on case level using information in %sample_info_file and %regexp
    case_qc(
        {
            evaluate_plink_gender => $evaluate_plink_gender,
            qc_data_href          => \%qc_data,
            qc_header_href        => \%qc_header,
            qc_recipe_data_href   => \%qc_recipe_data,
            regexp_href           => \%regexp,
            sample_info_href      => \%sample_info,
        }
    );

    evaluate_analysis(
        {
            eval_metric_file => $eval_metric_file,
            qc_data_href     => \%qc_data,
            sample_info_href => \%sample_info,
            skip_evaluation  => $skip_evaluation,
        }
    );

    ## Writes a qc data metrics to file
    store_metrics(
        {
            qc_data_href          => \%qc_data,
            sample_info_href      => \%sample_info,
            store_metrics_outfile => $store_metrics_outfile,
        }
    );

    parse_limit_qc_output(
        {
            limit_qc_output => $limit_qc_output,
            qc_href         => \%qc_data,
        }
    );

    ## Writes a qc data hash to file
    write_to_file(
        {
            data_href => \%qc_data,
            format    => q{yaml},
            path      => $outfile,
        }
    );
    $log->info( q{Wrote: } . $outfile );

    return;
}

######################
####Sub routines######
######################

sub case_qc {

## Function : Extracts all Qc data on case level using information in %sample_info_file and %regexp
## Returns  :
## Arguments: $evaluate_plink_gender => Evaluate plink gender
##          : $qc_data_href          => Qc data hash {REF}
##          : $qc_header_href        => Save header(s) in each outfile {REF}
##          : $qc_recipe_data_href   => Hash to save data in each outfile {REF}
##          : $regexp_href           => RegExp hash {REF}
##          : $sample_info_href      => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $evaluate_plink_gender;
    my $qc_data_href;
    my $qc_header_href;
    my $qc_recipe_data_href;
    my $regexp_href;
    my $sample_info_href;

    my $tmpl = {
        evaluate_plink_gender => {
            store       => \$evaluate_plink_gender,
            strict_type => 1,
        },
        qc_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_data_href,
            strict_type => 1,
        },
        qc_header_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_header_href,
            strict_type => 1,
        },
        qc_recipe_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_recipe_data_href,
            strict_type => 1,
        },
        regexp_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$regexp_href,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Qccollect qw{ plink_gender_check relation_check };
    use MIP::Qc_data
      qw{ parse_qc_recipe_data parse_regexp_hash_and_collect set_qc_data_recipe_info };
    use MIP::Sample_info qw{ get_sample_info_case_recipe_attributes };

  RECIPE:
    for my $recipe ( keys %{ $sample_info_href->{recipe} } ) {

        my %attribute = get_sample_info_case_recipe_attributes(
            {
                recipe_name      => $recipe,
                sample_info_href => $sample_info_href,
            }
        );

        my $outdirectory = $attribute{outdirectory};
        my $outfile      = $attribute{outfile};

        if ( exists $attribute{path} ) {

            ( $outfile, $outdirectory ) =
              fileparse( $attribute{path} );
        }

        ## Set package executable version from recipe to metrics hash
        set_qc_data_recipe_info(
            {
                key          => q{version},
                qc_data_href => $qc_data_href,
                recipe_name  => $recipe,
                value        => $attribute{version},
            }
        );

        ## Parses the RegExpHash structure to identify if the info is
        ## 1) Paragraf section(s) (both header and data line(s)
        ## 2) Seperate data line
        parse_regexp_hash_and_collect(
            {
                outdirectory        => $outdirectory,
                outfile             => $outfile,
                qc_recipe_data_href => $qc_recipe_data_href,
                qc_header_href      => $qc_header_href,
                recipe              => $recipe,
                regexp_href         => $regexp_href,
            }
        );

        ## Parse qc_recipe_data and extract information to qc_data
        parse_qc_recipe_data(
            {
                qc_data_href        => $qc_data_href,
                qc_header_href      => $qc_header_href,
                qc_recipe_data_href => $qc_recipe_data_href,
                recipe              => $recipe,
                regexp_href         => $regexp_href,
                sample_info_href    => $sample_info_href,
            }
        );

        ## Check gender for sample_id
        if (    $recipe eq q{plink_sexcheck}
            and $evaluate_plink_gender )
        {

            ## Check that assumed gender is supported by variants on chrX and chrY
            plink_gender_check(
                {
                    qc_data_href     => $qc_data_href,
                    sample_info_href => $sample_info_href,
                }
            );
        }
    }

    if (    exists $qc_data_href->{recipe}{pedigree_check}
        and exists $qc_data_href->{recipe}{relation_check} )
    {

        relation_check(
            {
                qc_data_href            => $qc_data_href,
                relationship_values_ref =>
                  \@{ $qc_data_href->{recipe}{relation_check}{sample_relation_check} },
                sample_info_href  => $sample_info_href,
                sample_orders_ref => \@{ $qc_data_href->{recipe}{pedigree_check}{sample_order} },
            }
        );
    }
    return;
}

sub sample_qc {

## Function : Collects all sample qc in files defined by sample_info_file and regular expressions defined by regexp.
## Returns  :
## Arguments: $qc_data_href        => Qc data hash {REF}
##          : $qc_header_href      => Save header(s) in each outfile {REF}
##          : $qc_recipe_data_href => Hash to save data in each outfile {REF}
##          : $regexp_href         => RegExp hash {REF}
##          : $sample_info_href    => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $qc_data_href;
    my $qc_header_href;
    my $qc_recipe_data_href;
    my $regexp_href;
    my $sample_info_href;

    my $tmpl = {
        qc_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_data_href,
            strict_type => 1,
        },
        qc_header_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_header_href,
            strict_type => 1,
        },
        qc_recipe_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_recipe_data_href,
            strict_type => 1,
        },
        regexp_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$regexp_href,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Qccollect qw{ chanjo_gender_check };
    use MIP::Qc_data qw{ parse_qc_recipe_data parse_regexp_hash_and_collect };
    use MIP::Sample_info qw{ get_sample_info_sample_recipe_attributes };

  SAMPLE_ID:
    for my $sample_id ( keys %{ $sample_info_href->{sample} } ) {

      RECIPE:
        for my $recipe ( keys %{ $sample_info_href->{sample}{$sample_id}{recipe} } ) {

          INFILE:
            for my $infile ( keys %{ $sample_info_href->{sample}{$sample_id}{recipe}{$recipe} } ) {

                my %attribute = get_sample_info_sample_recipe_attributes(
                    {
                        infile           => $infile,
                        recipe_name      => $recipe,
                        sample_id        => $sample_id,
                        sample_info_href => $sample_info_href,
                    }
                );

                my $outdirectory = $attribute{outdirectory};
                my $outfile      = $attribute{outfile};

                if ( exists $attribute{path} ) {

                    ( $outfile, $outdirectory ) =
                      fileparse( $attribute{path} );
                }

                ## Parses the RegExpHash structure to identify if the info is
                ## 1) Paragraf section(s) (both header and data line(s)
                ## 2) Seperate data line
                parse_regexp_hash_and_collect(
                    {
                        outdirectory        => $outdirectory,
                        outfile             => $outfile,
                        qc_header_href      => $qc_header_href,
                        qc_recipe_data_href => $qc_recipe_data_href,
                        recipe              => $recipe,
                        regexp_href         => $regexp_href,
                    }
                );

                ## Parse qc_recipe_data and extract information to qc_data
                parse_qc_recipe_data(
                    {
                        infile              => $infile,
                        qc_data_href        => $qc_data_href,
                        qc_header_href      => $qc_header_href,
                        qc_recipe_data_href => $qc_recipe_data_href,
                        recipe              => $recipe,
                        regexp_href         => $regexp_href,
                        sample_id           => $sample_id,
                        sample_info_href    => $sample_info_href,
                    }
                );

                ## Check gender for sample_id
                if ( $recipe eq q{chanjo_sexcheck} ) {

                    ## Check that assumed gender is supported by coverage on chrX and chrY
                    chanjo_gender_check(
                        {
                            infile           => $infile,
                            qc_data_href     => $qc_data_href,
                            recipe_name      => $recipe,
                            sample_id        => $sample_id,
                            sample_info_href => $sample_info_href,
                        }
                    );
                }
            }
        }
    }
    return;
}

sub parse_limit_qc_output {

## Function : Restrict output
## Returns  :
## Arguments: $limit_qc_output => Remove keys from regexp hash
##          : $qc_href         => Qccollect regexp hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $limit_qc_output;
    my $qc_href;

    my $tmpl = {
        limit_qc_output => {
            store       => \$limit_qc_output,
            strict_type => 1,
        },
        qc_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Delete surplus metrics
    delete $qc_href->{metrics};

    return if not $limit_qc_output;

    Readonly my @SKIP_QC_METRICS => qw{ variantevalall variantevalexome };

    foreach my $delete_metric_key (@SKIP_QC_METRICS) {

        _delete_key(
            {
                data_href  => $qc_href,
                delete_key => $delete_metric_key,
            }
        );
    }
    return;
}

sub _delete_key {

## Function : Delete key from nested hash
## Returns  :
## Arguments: $data_href  => Data {REF}
##          : $delete_key => Key to remove

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $data_href;
    my $delete_key;

    my $tmpl = {
        data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$data_href,
            strict_type => 1,
        },
        delete_key => {
            store       => \$delete_key,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  KEY_VALUE_PAIR:
    while ( my ( $key, $value ) = each %{$data_href} ) {

        if ( $key eq $delete_key ) {

            delete $data_href->{$delete_key};
            return;
        }

        elsif ( ref $data_href->{$key} eq q{HASH} ) {

            _delete_key(
                {
                    data_href  => $value,
                    delete_key => $delete_key,
                }
            );
        }
    }
    return;
}

1;
