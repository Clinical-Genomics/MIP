package MIP::Set::Pedigree;

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
use List::MoreUtils qw { any };
use Readonly;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ set_pedigree_capture_kit_info };
}

## Constants
Readonly my $SPACE => q{ };

sub set_pedigree_capture_kit_info {

## Function : Add and set capture kit for each individual
## Returns  :
## Arguments: $active_parameter_href   => Active parameters for this analysis hash {REF}
##          : $parameter_href          => Parameter hash {REF}
##          : $pedigree_href           => YAML pedigree info hash {REF}
##          : $sample_info_href        => Info on samples and case hash {REF}
##          : $user_supply_switch_href => The user supplied info switch {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $parameter_href;
    my $pedigree_href;
    my $sample_info_href;
    my $user_supply_switch_href;

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
        pedigree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$pedigree_href,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        user_supply_switch_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$user_supply_switch_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::Parameter qw{ get_capture_kit };

    my %exom_target_bed_file_tracker;

  SAMPLE_HREF:
    foreach my $pedigree_sample_href ( @{ $pedigree_href->{samples} } ) {

        ## Alias
        my $sample_id = $pedigree_sample_href->{sample_id};
        my $capture_kit =
          $sample_info_href->{sample}{$sample_id}{capture_kit};

        ## No recorded capture kit from pedigree or previous run
        next SAMPLE_HREF if ( not $capture_kit );

        ## Return a capture kit depending on user info
        my $exome_target_bed_file = get_capture_kit(
            {
                capture_kit => $capture_kit,
                supported_capture_kit_href =>
                  $parameter_href->{supported_capture_kit}{default},
                user_supplied_parameter_switch =>
                  $user_supply_switch_href->{exome_target_bed},
            }
        );

        ## No capture kit returned
        next SAMPLE_HREF if ( not $exome_target_bed_file );

        push @{ $exom_target_bed_file_tracker{$exome_target_bed_file} }, $sample_id;
    }

  BED_FILE:
    foreach my $exome_target_bed_file ( keys %exom_target_bed_file_tracker ) {

        ## We have read capture kits from pedigree and
        ## need to transfer to active_parameters
        $active_parameter_href->{exome_target_bed}{$exome_target_bed_file}
          = join q{,},
          @{ $exom_target_bed_file_tracker{$exome_target_bed_file} };
    }
    return;
}

1;
