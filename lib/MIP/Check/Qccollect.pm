package MIP::Check::Qccollect;

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
use MIP::Constants qw{ $COLON $NEWLINE $SPACE $TAB };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ chanjo_gender_check plink_gender_check relation_check };
}

sub chanjo_gender_check {

## Function : Checks that the gender predicted by chanjo sexcheck is confirmed in the pedigee for the sample
## Returns  :
## Arguments: $infile                 => Infile {REF}
##          : $qc_data_href           => Qc data hash {REF}
##          : $recipe_name            => Recipe to set attributes for
##          : $sample_id              => Sample ID
##          : $sample_info_href       => Info on samples and case hash {REF}

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

    use MIP::Get::Parameter qw{ get_pedigree_sample_id_attributes };
    use MIP::Qc_data qw{ get_qc_data_sample_recipe_attributes set_qc_data_recipe_info };

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
    use MIP::Get::Parameter qw{ get_pedigree_sample_id_attributes };

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

    use MIP::Qc_data qw{ set_qc_data_recipe_info };
    use MIP::Get::Parameter qw{ get_pedigree_sample_id_attributes };

    ## Constants
    Readonly my $RELATIONSHIP_CUTOFF => 0.70;

    ## Stores case relations and pairwise comparisons case{$sample_id}{$sample_id}["column"] -> [pairwise]
    my %case;
    my $incorrect_relation = 0;
    my $sample_id_counter  = 0;

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

    ## Father_id for the case
    my $father_id = q{YYY};

    ## Mother_id for the case
    my $mother_id = q{XXX};

    ## Collect father and mother id
  SAMPLE_ID:
    for my $sample_id ( keys %case ) {

        ## Currently only 1 father or mother per pedigree is supported

        ## Save father_id if not 0
        if ( $sample_info_href->{sample}{$sample_id}{father} ne 0 ) {

            $father_id = get_pedigree_sample_id_attributes(
                {
                    attribute        => q{father},
                    sample_id        => $sample_id,
                    sample_info_href => $sample_info_href,
                }
            );
        }

        ## Save mother_id if not 0
        if ( $sample_info_href->{sample}{$sample_id}{mother} ne 0 ) {

            $mother_id = get_pedigree_sample_id_attributes(
                {
                    attribute        => q{mother},
                    sample_id        => $sample_id,
                    sample_info_href => $sample_info_href,
                }
            );
        }
    }

  SAMPLE_ID:
    for my $sample_id ( keys %case ) {

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
