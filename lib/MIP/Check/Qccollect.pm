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
    our @EXPORT_OK = qw{ chanjo_gender_check plink_gender_check };
}

## Constants
Readonly my $FIELD_COUNTER => 2;

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

    use MIP::Get::Parameter qw{ get_pedigree_sample_id_attributes };

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

  SAMPLE_SEX:
    foreach
      my $data_metric ( @{ $qc_data_href->{recipe}{plink_sexcheck}{sample_sexcheck} } )
    {

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
            push @{ $qc_data_href->{recipe}{plink_gender_check} }, $sample_id . q{:PASS};
            next SAMPLE_SEX;
        }

        push @{ $qc_data_href->{recipe}{plink_gender_check} }, $sample_id . q{:FAIL};
    }
    return;
}

1;
