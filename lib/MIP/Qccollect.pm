package MIP::Qccollect;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use List::Util qw{ any };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $COLON $LOG_NAME $NEWLINE $SPACE $TAB $UNDERSCORE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      check_qc_metric
      chanjo_gender_check
      define_evaluate_metric
      evaluate_analysis
      evaluate_case_qc_parameters
      evaluate_sample_qc_parameters
      get_case_pairwise_comparison
      get_eval_expression
      get_parent_ids
      parse_sample_recipe_qc_metric
      parse_sample_qc_metric
      plink_gender_check
      plink_relation_check
      relation_check
      set_case_eval_metrics
      set_eval_expression
      set_sample_eval_metrics
      set_in_analysis_eval_metric
    };
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
        metric       => { defined => 1, required => 1, store => \$metric, strict_type => 1, },
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

    if (    defined $chanjo_sexcheck_gender
        and defined $sample_id_sex
        and exists $gender_map{$chanjo_sexcheck_gender}{$sample_id_sex} )
    {

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
      q{Status:} . $recipe_name . q{:} . $qc_data_href->{sample}{$sample_id}{$recipe_name};

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
## Arguments : eval_metric_file  => File with evaluation metrics
##           : $sample_info_href => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $eval_metric_file;
    my $sample_info_href;

    my $tmpl = {
        eval_metric_file => {
            defined     => 1,
            required    => 1,
            store       => \$eval_metric_file,
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

    use Clone qw{ clone };
    use MIP::Io::Read qw{ read_from_file };

    my %eval_metric = read_from_file(
        {
            format => q{yaml},
            path   => $eval_metric_file,
        }
    );

    ## Hash to store relevant metrics
    my %analysis_eval_metric;

    ## Case level
    set_case_eval_metrics(
        {
            analysis_eval_metric_href => \%analysis_eval_metric,
            eval_metric_href          => \%eval_metric,
            sample_info_href          => $sample_info_href,
        }
    );

    ## Sample level
  SAMPLE_ID:
    foreach my $sample_id ( keys %{ $sample_info_href->{sample} } ) {

        my %clone_eval_metric = %{ clone( \%eval_metric ) };

        set_sample_eval_metrics(
            {
                analysis_eval_metric_href => \%analysis_eval_metric,
                eval_metric_href          => \%clone_eval_metric,
                sample_id                 => $sample_id,
                sample_info_href          => $sample_info_href,
            }
        );
    }
    return %analysis_eval_metric;
}

sub evaluate_analysis {

## Function  : Read in evaluation metrics and evaluate on sample and case level
## Returns   :
## Arguments : eval_metric_file  => File with evaluation metrics
##           : qc_data_href      => QC data {REF}
##           : $sample_info_href => Info on samples and case hash {REF}
##           : $skip_evaluation  => Skip evaluation

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $eval_metric_file;
    my $qc_data_href;
    my $sample_info_href;
    my $skip_evaluation;

    my $tmpl = {
        eval_metric_file => {
            store       => \$eval_metric_file,
            strict_type => 1,
        },
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
        skip_evaluation => {
            allow       => [ undef, 0, 1 ],
            required    => 1,
            store       => \$skip_evaluation,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    if ($skip_evaluation) {

        $log->info(q{Skipping QC evaluation});
        return;
    }

    ## Defines recipes, metrics and thresholds to evaluate
    my %evaluate_metric = define_evaluate_metric(
        {
            eval_metric_file => $eval_metric_file,
            sample_info_href => $sample_info_href,
        }
    );

    ## Exit if there aren't any evaluation metrics
    if ( not %evaluate_metric ) {
        $log->warn(
qq{Supplied eval_metric_file, $eval_metric_file, contains no evaluation metrics that are relevant to the analysis}
        );
        $log->warn(qq{Skipping QC evaluation});
        return;
    }

    ## Evaluate the metrics
    evaluate_case_qc_parameters(
        {
            evaluate_metric_href => \%evaluate_metric,
            qc_data_href         => $qc_data_href,
        }
    );

    evaluate_sample_qc_parameters(
        {
            evaluate_metric_href => \%evaluate_metric,
            qc_data_href         => $qc_data_href,
        }
    );

    return 1;
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
        my @pairwise_comparisons = splice @relationship_values, 0, scalar @{$sample_orders_ref};

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

        ## Return father_id and mother_id if defined and not 0
        next SAMPLE_ID if ( not defined $father_id and not defined $mother_id );
        next SAMPLE_ID if ( $father_id eq 0        and $mother_id eq 0 );

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
                    metric                => $metric,
                    qc_data_href          => $qc_data_href,
                    qc_metric_value       => $qc_data_recipe_href->{$metric},
                    recipe                => $recipe_name,
                    reference_metric_href =>
                      $evaluate_metric_href->{$sample_id}{$recipe_name}{$metric},
                }
            );
            next METRIC;
        }

        ## No header data for metric
        next METRIC if ( not exists $qc_data_recipe_href->{header} );

      HEADER:
        for my $data_header ( keys %{ $qc_data_recipe_href->{header} } ) {

            ## Metric does not exist in header
            next HEADER
              if ( not exists $qc_data_recipe_href->{header}{$data_header}{$metric} );

            check_qc_metric(
                {
                    metric                => $metric,
                    qc_data_href          => $qc_data_href,
                    qc_metric_value       => $qc_data_recipe_href->{header}{$data_header}{$metric},
                    recipe                => $recipe_name,
                    reference_metric_href =>
                      $evaluate_metric_href->{$sample_id}{$recipe_name}{$metric},
                }
            );
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

        if (    defined $plink_sexcheck_gender
            and defined $sample_id_sex
            and exists $gender_map{$plink_sexcheck_gender}{$sample_id_sex} )
        {

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

## Function : Uses the .mibs file produced by PLINK to test if case members are indeed related
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
            required    => 1,
            store       => \$sample_orders_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Qccollect qw{ get_parent_ids };
    use MIP::Qc_data qw{ set_qc_data_recipe_info };
    use MIP::Sample_info qw{ get_pedigree_sample_id_attributes };

    return if ( not @{$relationship_values_ref} and not @{$sample_orders_ref} );

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
                if ( $relative_metric eq 1 ) {

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
    return 1;
}

sub set_case_eval_metrics {

## Function : Set evaluation metrics on case level
## Returns  :
## Arguments: $analysis_eval_metric_href => Hash with evaluation metrics for the analysis {ref}
##          : $eval_metric_href          => Hash with evaluation metrics {ref}
##          : $sample_info_href          => Hash with sample info {ref}

    my ($arg_href) = @_;

    ## flatten argument(s)
    my $analysis_eval_metric_href;
    my $eval_metric_href;
    my $sample_info_href;

    my $tmpl = {
        analysis_eval_metric_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$analysis_eval_metric_href,
            strict_type => 1,
        },
        eval_metric_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$eval_metric_href,
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

    check( $tmpl, $arg_href, 1 ) or croak q{could not parse arguments!};

    use MIP::Analysis qw{ get_overall_analysis_type };

    my $consensus_analysis_type = get_overall_analysis_type(
        {
            analysis_type_href => $sample_info_href->{analysis_type},
        }
    );
    my $pipeline_eval_metric_href = $eval_metric_href->{$consensus_analysis_type};

  RECIPE_OUTFILE:
    foreach my $recipe ( keys %{ $sample_info_href->{recipe} } ) {

        if ( any { $_ eq $recipe } ( keys %{$pipeline_eval_metric_href} ) ) {

            my %eval_expression = get_eval_expression(
                {
                    eval_metric_href => $pipeline_eval_metric_href,
                    recipe           => $recipe,
                }
            );

            set_eval_expression(
                {
                    analysis_eval_metric_href => $analysis_eval_metric_href,
                    recipe                    => $recipe,
                    eval_expression_href      => \%eval_expression,
                }
            );
        }
    }
    return;
}

sub set_sample_eval_metrics {

## function : Set evaluation metrics on sample level
## returns  :
## arguments: $analysis_eval_metric_href => Hash with evaluation metrics for the analysis {ref}
##          : $eval_metric_href          => Hash with evaluation metrics {ref}
##          : $sample_id                 => Sample id
##          : $sample_info_href          => Hash with sample info {ref}

    my ($arg_href) = @_;

    ## flatten argument(s)
    my $analysis_eval_metric_href;
    my $eval_metric_href;
    my $sample_id;
    my $sample_info_href;

    my $tmpl = {
        analysis_eval_metric_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$analysis_eval_metric_href,
            strict_type => 1,
        },
        eval_metric_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$eval_metric_href,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{could not parse arguments!};

    use MIP::Sample_info qw{ get_pedigree_sample_id_attributes };

    ## Get info from sample info hash
    my %sample_attribute = get_pedigree_sample_id_attributes(
        {
            sample_id        => $sample_id,
            sample_info_href => $sample_info_href,
        }
    );

    my $analysis_type             = $sample_attribute{analysis_type};
    my $expected_coverage         = $sample_attribute{expected_coverage};
    my %recipe                    = %{ $sample_attribute{recipe} };
    my $pipeline_eval_metric_href = $eval_metric_href->{$analysis_type};

  RECIPE_OUTFILE:
    foreach my $recipe ( keys %recipe ) {

        if ( any { $_ eq $recipe } ( keys %{$pipeline_eval_metric_href} ) ) {

            my %eval_expression = get_eval_expression(
                {
                    eval_metric_href => $pipeline_eval_metric_href,
                    recipe           => $recipe,
                }
            );

            set_eval_expression(
                {
                    analysis_eval_metric_href => $analysis_eval_metric_href,
                    recipe                    => $recipe,
                    eval_expression_href      => \%eval_expression,
                    sample_id                 => $sample_id,
                }
            );
        }
    }

    ## Special case foir expected coverage
    if ($expected_coverage) {

        my %eval_expression = (
            MEDIAN_TARGET_COVERAGE => {
                lt => $expected_coverage,
            },
        );
        set_eval_expression(
            {
                analysis_eval_metric_href => $analysis_eval_metric_href,
                recipe                    => q{collecthsmetrics},
                eval_expression_href      => \%eval_expression,
                sample_id                 => $sample_id,
            }
        );
    }

    return;
}

sub set_eval_expression {

## function : Set evaluation expression in hash
## returns  :
## arguments: $analysis_eval_metric_href => Hash with evaluation metrics for the analysis {ref}
##          : $eval_expression_href      => Hash with evaluation expression {ref}
##          : $recipe                    => Recipe
##          : $sample_id                 => Sample id

    my ($arg_href) = @_;

    ## flatten argument(s)
    my $analysis_eval_metric_href;
    my $eval_expression_href;
    my $recipe;
    my $sample_id;

    my $tmpl = {
        analysis_eval_metric_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$analysis_eval_metric_href,
            strict_type => 1,
        },
        eval_expression_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$eval_expression_href,
            strict_type => 1,
        },
        recipe => {
            defined     => 1,
            required    => 1,
            store       => \$recipe,
            strict_type => 1,
        },
        sample_id => {
            store       => \$sample_id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{could not parse arguments!};

    if ($sample_id) {

        $analysis_eval_metric_href->{$sample_id}{$recipe} = $eval_expression_href;
        return;
    }

    $analysis_eval_metric_href->{$recipe} = $eval_expression_href;
    return;
}

sub get_eval_expression {

## function : Get evaluation expression
## returns  : %recipe_eval_metric
## arguments: $eval_metric_href => Hash with evaluation metrics for the analysis {ref}
##          : $recipe           => Recipe
##          : $sample_id        => Sample id

    my ($arg_href) = @_;

    ## flatten argument(s)
    my $eval_metric_href;
    my $recipe;
    my $sample_id;

    my $tmpl = {
        eval_metric_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$eval_metric_href,
            strict_type => 1,
        },
        recipe => {
            defined     => 1,
            required    => 1,
            store       => \$recipe,
            strict_type => 1,
        },
        sample_id => {
            store       => \$sample_id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{could not parse arguments!};

    if ($sample_id) {

        return %{ $eval_metric_href->{$sample_id}{$recipe} };
    }

    return %{ $eval_metric_href->{$recipe} };
}

1;
