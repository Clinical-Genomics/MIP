package MIP::Qccollect;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $COLON $NEWLINE $SPACE $TAB $UNDERSCORE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      check_qc_metric
      chanjo_gender_check
      define_evaluate_metric
      evaluate_case_qc_parameters
      evaluate_sample_qc_parameters
      get_case_pairwise_comparison
      get_parent_ids
      parse_sample_recipe_qc_metric
      parse_sample_qc_metric
      plink_gender_check
      plink_relation_check
      relation_check };
}

sub check_qc_metric {

## Function : Check and add result of check to qc data hash
## Returns  :
## Arguments: $metric                => Metric to evaluate
##          : $qc_data_href          => Qc data hash {REF}
##          : $qc_metric_value       => Qc metric value
##          : $recipe                => The recipe to examine
##          : $reference_metric_href => Metrics to evaluate

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $metric;
    my $qc_data_href;
    my $qc_metric_value;
    my $recipe;
    my $reference_metric_href;

    my $tmpl = {
        metric => { defined => 1, required => 1, store => \$metric, strict_type => 1, },
        qc_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_data_href,
            strict_type => 1,
        },
        qc_metric_value => {
            defined     => 1,
            required    => 1,
            store       => \$qc_metric_value,
            strict_type => 1,
        },
        recipe => { defined => 1, required => 1, store => \$recipe, strict_type => 1, },
        reference_metric_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$reference_metric_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Qc_data qw{ add_qc_data_evaluation_info };

    my $status = q{FAILED:};

    if ( exists $reference_metric_href->{lt} ) {

        ## Determine status - if lower than add to hash. otherwise PASS and do not include
        if ( $qc_metric_value < $reference_metric_href->{lt} ) {

            ## Add to status string and then to hash
            $status .= $recipe . $UNDERSCORE . $metric . $COLON . $qc_metric_value;
            add_qc_data_evaluation_info(
                {
                    qc_data_href => $qc_data_href,
                    recipe_name  => $recipe,
                    value        => $status,
                }
            );
        }
    }

    if ( exists $reference_metric_href->{gt} ) {

        ## Determine status - if greater than add to hash. otherwise PASS and do not include
        if ( $qc_metric_value > $reference_metric_href->{gt} ) {

            $status .= $recipe . $UNDERSCORE . $metric . $COLON . $qc_metric_value;
            add_qc_data_evaluation_info(
                {
                    qc_data_href => $qc_data_href,
                    recipe_name  => $recipe,
                    value        => $status,
                }
            );
        }
    }
    return;
}

sub chanjo_gender_check {

## Function : Checks that the gender predicted by chanjo sexcheck is confirmed in the pedigee for the sample
## Returns  :
## Arguments: $infile           => Infile {REF}
##          : $qc_data_href     => Qc data hash {REF}
##          : $recipe_name      => Recipe to set attributes for
##          : $sample_id        => Sample ID
##          : $sample_info_href => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile;
    my $qc_data_href;
    my $recipe_name;
    my $sample_id;
    my $sample_info_href;

    my $tmpl = {
        infile => {
            defined     => 1,
            required    => 1,
            store       => \$infile,
            strict_type => 1,
        },
        qc_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_data_href,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
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

    use MIP::Qc_data qw{ get_qc_data_sample_recipe_attributes set_qc_data_recipe_info };
    use MIP::Sample_info qw{ get_pedigree_sample_id_attributes };

    ## Get chanjo estimated gender
    my $chanjo_sexcheck_gender = get_qc_data_sample_recipe_attributes(
        {
            attribute    => q{gender},
            infile       => $infile,
            recipe_name  => $recipe_name,
            sample_id    => $sample_id,
            qc_data_href => $qc_data_href,
        }
    );

    ## Get sample id sex
    my $sample_id_sex = get_pedigree_sample_id_attributes(
        {
            attribute        => q{sex},
            sample_id        => $sample_id,
            sample_info_href => $sample_info_href,
        }
    );

    ## Create map of allowed keys per sex
    my %gender_map = (
        female => {
            female => undef,
            2      => undef,
        },
        male => {
            male => undef,
            1    => undef,
        },
        other => {
            other => undef,
            0     => undef,
        },
        unknown => {
            unknown => undef,
            0       => undef,
        },
    );

## Set initial gender check status
    my $status = q{FAIL};

    if ( exists $gender_map{$chanjo_sexcheck_gender}{$sample_id_sex} ) {

        $status = q{PASS};
    }
    set_qc_data_recipe_info(
        {
            key          => q{gender_check},
            infile       => $infile,
            qc_data_href => $qc_data_href,
            recipe_name  => $recipe_name,
            sample_id    => $sample_id,
            value        => $status,
        }
    );
    return;
}

sub evaluate_sample_qc_parameters {

## Function : Evaluate sample qc metrics to detect metrics falling below threshold
## Returns  :
## Arguments: $evaluate_metric_href => Hash for metrics to evaluate
##          : $qc_data_href         => QC data hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $evaluate_metric_href;
    my $qc_data_href;

    my $tmpl = {
        evaluate_metric_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$evaluate_metric_href,
            strict_type => 1,
        },
        qc_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_data_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Qccollect qw{ parse_sample_qc_metric plink_relation_check };

  SAMPLE_LEVEL:
    for my $sample_id ( keys %{ $qc_data_href->{sample} } ) {

      INFILE:
        for my $infile ( keys %{ $qc_data_href->{sample}{$sample_id} } ) {

            ## Skip evaluation for these infiles reg exp matches
            next INFILE if ( $infile =~ / Undetermined | evaluation /ixsm );

            ## Need to check and possibly add status for plink relation check
            if ( $infile =~ /relation_check/isxm ) {

                plink_relation_check(
                    {
                        qc_data_href => $qc_data_href,
                        recipe_name  => $infile,
                        sample_id    => $sample_id,
                    }
                );
                next INFILE;
            }

            parse_sample_recipe_qc_metric(
                {
                    evaluate_metric_href => $evaluate_metric_href,
                    qc_data_href         => $qc_data_href,
                    infile               => $infile,
                    sample_id            => $sample_id,
                }
            );

        }
    }
    return 1;
}

sub plink_relation_check {

## Function : Evaluate plink relation check not equal PASS
## Returns  :
## Arguments: $qc_data_href => QC data hash {REF}
##          : $recipe_name  => Plink recipe name for relation check
##          : $sample_id    => Sample ID

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $qc_data_href;
    my $recipe_name;
    my $sample_id;

    my $tmpl = {
        qc_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_data_href,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Qc_data qw{ add_qc_data_evaluation_info };

    return if ( not exists $qc_data_href->{sample}{$sample_id}{$recipe_name} );

    return if ( $qc_data_href->{sample}{$sample_id}{$recipe_name} eq q{PASS} );

    ## Set FAIL status
    my $status =
        q{Status:}
      . $recipe_name . q{:}
      . $qc_data_href->{sample}{$sample_id}{$recipe_name};

    ## Add to QC data at case level
    add_qc_data_evaluation_info(
        {
            qc_data_href => $qc_data_href,
            recipe_name  => $recipe_name,
            value        => $status,
        }
    );
    return;
}

sub define_evaluate_metric {

## Function  : Sets recipes, metrics and thresholds to be evaluated
## Returns   :
## Arguments : $sample_info_href => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_info_href;

    my $tmpl = {
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Sample_info qw{ get_pedigree_sample_id_attributes };

    ## Constants
    Readonly my $FRACTION_DUPLICATES         => 0.2;
    Readonly my $FRACTION_OF_COMMON_VARIANTS => 0.55;
    Readonly my $FRACTION_OF_ERRORS          => 0.06;
    Readonly my $PCT_ADAPTER                 => 0.0005;
    Readonly my $PCT_PF_READS_ALIGNED        => 0.95;
    Readonly my $PCT_TARGET_BASES_10X        => 0.95;
    Readonly my $PERCENTAGE_MAPPED_READS     => 95;

    my %evaluate_metric;

  SAMPLE_ID:
    foreach my $sample_id ( keys %{ $sample_info_href->{sample} } ) {

        $evaluate_metric{$sample_id}{bamstats}{percentage_mapped_reads}{lt} =
          $PERCENTAGE_MAPPED_READS;
        $evaluate_metric{$sample_id}{collecthsmetrics}{PCT_TARGET_BASES_10X}{lt} =
          $PCT_TARGET_BASES_10X;
        $evaluate_metric{$sample_id}{collectmultiplemetrics}{PCT_PF_READS_ALIGNED}{lt} =
          $PCT_PF_READS_ALIGNED;
        $evaluate_metric{$sample_id}{collectmultiplemetrics}{PCT_ADAPTER}{gt} =
          $PCT_ADAPTER;
        $evaluate_metric{$sample_id}{markduplicates}{fraction_duplicates}{gt} =
          $FRACTION_DUPLICATES;
        $evaluate_metric{variant_integrity_ar_mendel}{fraction_of_errors}{gt} =
          $FRACTION_OF_ERRORS;
        $evaluate_metric{variant_integrity_ar_father}{fraction_of_common_variants}{lt} =
          $FRACTION_OF_COMMON_VARIANTS;

        ## Get sample id expected_coverage
        my $expected_coverage = get_pedigree_sample_id_attributes(
            {
                attribute        => q{expected_coverage},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );
        if ($expected_coverage) {

            $evaluate_metric{$sample_id}{collecthsmetrics}{MEAN_TARGET_COVERAGE}{lt} =
              $expected_coverage;
        }
    }
    return %evaluate_metric;
}

sub evaluate_case_qc_parameters {

## Function : Evaluate case qc metrics to detect metrics falling below threshold
## Returns  :
## Arguments: $evaluate_metric_href => Hash for metrics to evaluate
##          : $qc_data_href         => QC data hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $evaluate_metric_href;
    my $qc_data_href;

    my $tmpl = {
        qc_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_data_href,
            strict_type => 1,
        },
        evaluate_metric_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$evaluate_metric_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Qccollect qw{ check_qc_metric };

  RECIPE:
    for my $recipe ( keys %{$evaluate_metric_href} ) {

        next RECIPE if ( not exists $qc_data_href->{recipe}{$recipe} );

      METRIC:
        for my $metric ( keys %{ $evaluate_metric_href->{$recipe} } ) {

            next METRIC if ( not exists $qc_data_href->{recipe}{$recipe}{$metric} );

            check_qc_metric(
                {
                    metric                => $metric,
                    qc_data_href          => $qc_data_href,
                    qc_metric_value       => $qc_data_href->{recipe}{$recipe}{$metric},
                    recipe                => $recipe,
                    reference_metric_href => $evaluate_metric_href->{$recipe}{$metric},
                }
            );
        }
    }
    return 1;
}

sub get_case_pairwise_comparison {

## Function : Uses the .mibs file produced by PLINK to test if case members are indeed related.
## Returns  : %case
## Arguments: $relationship_values_ref => All relationship estimations {REF}
##          : $sample_orders_ref       => The sample order so that correct estimation can be connected to the correct sample_ids {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $relationship_values_ref;
    my $sample_orders_ref;

    my $tmpl = {
        relationship_values_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$relationship_values_ref,
            strict_type => 1,
        },
        sample_orders_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$sample_orders_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %case;
    my $sample_id_counter = 0;

    ## Copy array to avoid removing actual values in downstream splice
    my @relationship_values = @{$relationship_values_ref};

    ## Splice all relationship estimations from regexp into pairwise comparisons calculated for each sample_id
  RELATIONSHIP:
    foreach my $relationship (@relationship_values) {

        ## Splices array into each sample_ids line
        my @pairwise_comparisons = splice @relationship_values, 0,
          scalar @{$sample_orders_ref};

        ## All columns in .mibs file
      COLUMN:
        while ( my ( $column_index, $sample_to_compare ) = each @{$sample_orders_ref} ) {

            ## Store sample_id, case membersID (including self) and each pairwise comparison.
            ## Uses array to accomodate sibling info.
            my $sample = $sample_orders_ref->[$sample_id_counter];
            push
              @{ $case{$sample}{$sample_to_compare} },
              $pairwise_comparisons[$column_index];
        }
        ## Increment counter for next sample to use as base in comparisons
        $sample_id_counter++;
    }
    return %case;
}

sub get_parent_ids {

## Function : Get parent ids for case and sample info hash
## Returns  : father_id, mother_id
## Arguments: $case_href        => Case relations and comparisons
##          : $sample_info_href => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $case_href;
    my $sample_info_href;

    my $tmpl = {
        case_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$case_href,
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

    use MIP::Sample_info qw{ get_pedigree_sample_id_attributes };

  SAMPLE_ID:
    for my $sample_id ( keys %{$case_href} ) {

        ## Currently only 1 father or mother per pedigree is supported

        my $father_id = get_pedigree_sample_id_attributes(
            {
                attribute        => q{father},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );
        my $mother_id = get_pedigree_sample_id_attributes(
            {
                attribute        => q{mother},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        ## Save father_id and  mother_id if not 0
        next SAMPLE_ID if ( $father_id eq 0 and $mother_id eq 0 );

        return $father_id, $mother_id;
    }
    ## Return fake ids to allow sibling comparisons for pedigrees without parents
    return q{YYY_father}, q{XXX_mother};
}

sub parse_sample_recipe_qc_metric {

## Function : Parse sample recipe qc data to check metrics
## Returns  :
## Arguments: $evaluate_metric_href => Hash for metrics to evaluate
##          : $infile               => Infile name
##          : $qc_data_href         => QC data hash {REF}
##          : $sample_id            => Sample ID

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $evaluate_metric_href;
    my $infile;
    my $qc_data_href;
    my $sample_id;

    my $tmpl = {
        evaluate_metric_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$evaluate_metric_href,
            strict_type => 1,
        },
        infile => {
            defined     => 1,
            required    => 1,
            store       => \$infile,
            strict_type => 1,
        },
        qc_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_data_href,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  RECIPE:
    for my $recipe ( keys %{ $qc_data_href->{sample}{$sample_id}{$infile} } ) {

        next RECIPE
          if ( not exists $evaluate_metric_href->{$sample_id}{$recipe} );

        ## Alias qc data for recipe
        my $qc_data_recipe_href =
          \%{ $qc_data_href->{sample}{$sample_id}{$infile}{$recipe} };

        ## Parse sample recipe qc data and/or header to check metrics
        parse_sample_qc_metric(
            {
                evaluate_metric_href => $evaluate_metric_href,
                qc_data_href         => $qc_data_href,
                qc_data_recipe_href  => $qc_data_recipe_href,
                recipe_name          => $recipe,
                sample_id            => $sample_id,
            }
        );

    }
    return 1;
}

sub parse_sample_qc_metric {

## Function : Parse sample recipe qc data and/or header to check metrics
## Returns  :
## Arguments: $evaluate_metric_href => Hash for metrics to evaluate
##          : $qc_data_href         => Qc data hash {REF}
##          : $qc_data_recipe_href  => Recipe specific qc data
##          : $recipe_name          => Recipe name
##          : $sample_id            => Sample ID

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $evaluate_metric_href;
    my $qc_data_href;
    my $qc_data_recipe_href;
    my $recipe_name;
    my $sample_id;

    my $tmpl = {
        evaluate_metric_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$evaluate_metric_href,
            strict_type => 1,
        },
        qc_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_data_href,
            strict_type => 1,
        },
        qc_data_recipe_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_data_recipe_href,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  METRIC:
    for my $metric ( keys %{ $evaluate_metric_href->{$sample_id}{$recipe_name} } ) {

        ## If metrics exists in qc data recipe
        if ( exists $qc_data_recipe_href->{$metric} ) {

            check_qc_metric(
                {
                    metric          => $metric,
                    qc_data_href    => $qc_data_href,
                    qc_metric_value => $qc_data_recipe_href->{$metric},
                    recipe          => $recipe_name,
                    reference_metric_href =>
                      $evaluate_metric_href->{$sample_id}{$recipe_name}{$metric},
                }
            );
            next METRIC;
        }

        ## No header data for metric
        next METRIC
          if ( not exists $qc_data_recipe_href->{header} );

      HEADER:
        for my $data_header ( keys %{ $qc_data_recipe_href->{header} } ) {

            ## Metric does not exist in header
            next HEADER
              if ( not exists $qc_data_recipe_href->{header}{$data_header}{$metric} );

            check_qc_metric(
                {
                    metric       => $metric,
                    qc_data_href => $qc_data_href,
                    qc_metric_value =>
                      $qc_data_recipe_href->{header}{$data_header}{$metric},
                    recipe => $recipe_name,
                    reference_metric_href =>
                      $evaluate_metric_href->{$sample_id}{$recipe_name}{$metric},
                }
            );
            next METRIC;
        }
    }
    return 1;
}

sub plink_gender_check {

## Function : Checks that the gender predicted by Plink sexcheck is confirmed in the pedigee for the sample
## Returns  :
## Arguments: $qc_data_href     => Qc data hash {REF}
##          : $sample_info_href => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $qc_data_href;
    my $sample_info_href;

    my $tmpl = {
        qc_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_data_href,
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

    use MIP::Qc_data qw{ add_qc_data_recipe_info get_qc_data_case_recipe_attributes };
    use MIP::Sample_info qw{ get_pedigree_sample_id_attributes };

    ## Constants
    Readonly my $FIELD_COUNTER => 2;

    ## Create map of allowed keys per sex
    my %gender_map = (
        2 => {
            female => undef,
            2      => undef,
        },
        1 => {
            male => undef,
            1    => undef,
        },
        0 => {
            other   => undef,
            0       => undef,
            unknown => undef,
        },
    );

    my $data_metrics_ref = get_qc_data_case_recipe_attributes(
        {
            attribute    => q{sample_sexcheck},
            qc_data_href => \%{$qc_data_href},
            recipe_name  => q{plink_sexcheck},
        }
    );

  SAMPLE_SEX:
    foreach my $data_metric ( @{$data_metrics_ref} ) {

        ## Array
        my ( $sample_id, $plink_sexcheck_gender, $unexpected_data ) = split $COLON,
          $data_metric, $FIELD_COUNTER + 1;

        ## Make sure that we get what we expect
        if ( defined $unexpected_data ) {

            carp q{Unexpected trailing garbage in data metric '} . $data_metric . q{':},
              $NEWLINE . $TAB . $unexpected_data . $NEWLINE;
        }

        ## Get sample id sex
        my $sample_id_sex = get_pedigree_sample_id_attributes(
            {
                attribute        => q{sex},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        if ( exists $gender_map{$plink_sexcheck_gender}{$sample_id_sex} ) {

            add_qc_data_recipe_info(
                {
                    key          => q{plink_gender_check},
                    qc_data_href => $qc_data_href,
                    recipe_name  => q{plink_sexcheck},
                    value        => $sample_id . q{:PASS},
                }
            );
            next SAMPLE_SEX;
        }

        add_qc_data_recipe_info(
            {
                key          => q{plink_gender_check},
                qc_data_href => $qc_data_href,
                recipe_name  => q{plink_sexcheck},
                value        => $sample_id . q{:FAIL},
            }
        );
    }
    return;
}

sub relation_check {

## Function : Uses the .mibs file produced by PLINK to test if case members are indeed related.
## Returns  :
## Arguments: $qc_data_href            => Qc data hash {REF}
##          : $relationship_values_ref => All relationship estimations {REF}
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $sample_orders_ref       => The sample order so that correct estimation can be connected to the correct sample_ids {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $qc_data_href;
    my $relationship_values_ref;
    my $sample_info_href;
    my $sample_orders_ref;

    my $tmpl = {
        qc_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qc_data_href,
            strict_type => 1,
        },
        relationship_values_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$relationship_values_ref,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        sample_orders_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$sample_orders_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Qccollect qw{ get_parent_ids };
    use MIP::Qc_data qw{ set_qc_data_recipe_info };
    use MIP::Sample_info qw{ get_pedigree_sample_id_attributes };

    ## Constants
    Readonly my $RELATIONSHIP_CUTOFF => 0.70;

    ## Stores case relations and pairwise comparisons case{$sample_id}{$sample_id}["column"] -> [pairwise]
    my %case = get_case_pairwise_comparison(
        {
            relationship_values_ref => $relationship_values_ref,
            sample_orders_ref       => $sample_orders_ref,
        }
    );

    my ( $father_id, $mother_id ) = get_parent_ids(
        {
            case_href        => \%case,
            sample_info_href => $sample_info_href,
        }
    );

  SAMPLE_ID:
    for my $sample_id ( keys %case ) {

        ## For each base sample
        my $incorrect_relation = 0;

      MEMBER:
        for my $members ( keys %{ $case{$sample_id} } ) {
            ## For every relation within case (mother/father/child)

          RELATIVE:
            foreach my $relative_metric ( @{ $case{$sample_id}{$members} } ) {

                ## Should only hit self
                if ( $relative_metric == 1 ) {

                    ## If self
                    next RELATIVE if ( $sample_id eq $members );

                    $incorrect_relation++;
                    set_qc_data_recipe_info(
                        {
                            key          => q{relation_check},
                            qc_data_href => $qc_data_href,
                            sample_id    => $sample_id,
                            value        => q{FAIL: Duplicated sample?},
                        }
                    );
                }
                elsif ( $relative_metric >= $RELATIONSHIP_CUTOFF ) {
                    ## Should include parent to child and child to siblings unless inbreed parents

                    ## Check comparison is not done between parents
                    next RELATIVE
                      if ( $sample_id ne $father_id && $sample_id ne $mother_id );

                    next RELATIVE if ( $members ne $father_id && $members ne $mother_id );

                    $incorrect_relation++;
                    set_qc_data_recipe_info(
                        {
                            key          => q{relation_check},
                            qc_data_href => $qc_data_href,
                            sample_id    => $sample_id,
                            value        => q{FAIL: Parents related?},
                        }
                    );
                }
                elsif ( $relative_metric < $RELATIONSHIP_CUTOFF ) {
                    ## Parents unless inbreed

                    next RELATIVE
                      if (  $sample_id eq $father_id
                        and $members eq $mother_id );

                    next RELATIVE
                      if (  $sample_id eq $mother_id
                        and $members eq $father_id );

                    $incorrect_relation++;
                    set_qc_data_recipe_info(
                        {
                            key          => q{relation_check},
                            qc_data_href => $qc_data_href,
                            sample_id    => $sample_id,
                            value        => qq{FAIL:$sample_id not related to $members},
                        }
                    );
                }
            }
        }
        if ( not $incorrect_relation ) {

            set_qc_data_recipe_info(
                {
                    key          => q{relation_check},
                    qc_data_href => $qc_data_href,
                    sample_id    => $sample_id,
                    value        => q{PASS},
                }
            );
        }
    }
    return;
}

1;
