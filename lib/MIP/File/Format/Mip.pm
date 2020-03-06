package MIP::File::Format::Mip;

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

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ build_file_prefix_tag fastq_file_name_regexp };
}

sub build_file_prefix_tag {

## Function : Build the file tags depending on which modules are used by the user to relevant chain.
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $case_id             => Family id {REF}
##          : $file_info_href        => Info on files hash {REF}
##          : $order_recipes_ref     => Order of addition to parameter array {REF}
##          : $parameter_href        => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $order_recipes_ref;
    my $parameter_href;

    ## Default(s)
    my $case_id;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        case_id => {
            default     => $arg_href->{active_parameter_href}{case_id},
            store       => \$case_id,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        order_recipes_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$order_recipes_ref,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Io::Recipes qw{ set_file_prefix_tag };

    ## Used to enable seqential build-up of file_tags between analysis recipes
    my %temp_file_ending;

  PARAMETER:
    foreach my $recipe_name ( @{$order_recipes_ref} ) {

        ## Alias
        my $current_chain = $parameter_href->{$recipe_name}{chain};
        my $file_tag      = $parameter_href->{$recipe_name}{file_tag};

        ## Only active parameters
        next PARAMETER
          if ( not defined $active_parameter_href->{$recipe_name} );

        ## Skip parameters with no file tag
        next PARAMETER if ( $file_tag eq q{nofile_tag} );

        ### Per sample_id
      SAMPLE_ID:
        foreach my $sample_id ( @{ $active_parameter_href->{sample_ids} } ) {

            ## Set file tag for recipe and update sequential build-up of fileending
            $temp_file_ending{$current_chain}{$sample_id} = set_file_prefix_tag(
                {
                    current_chain         => $current_chain,
                    file_tag              => $file_tag,
                    file_info_href        => $file_info_href,
                    id                    => $sample_id,
                    is_active_recipe      => $active_parameter_href->{$recipe_name},
                    recipe_name           => $recipe_name,
                    temp_file_ending_href => \%temp_file_ending,
                }
            );
        }

        ### Per case_id
        ## Set file tag for recipe and update sequential build-up of fileending
        $temp_file_ending{$current_chain}{$case_id} = set_file_prefix_tag(
            {
                current_chain         => $current_chain,
                file_tag              => $file_tag,
                file_info_href        => $file_info_href,
                id                    => $case_id,
                is_active_recipe      => $active_parameter_href->{$recipe_name},
                recipe_name           => $recipe_name,
                temp_file_ending_href => \%temp_file_ending,
            }
        );
    }
    return;
}

sub fastq_file_name_regexp {

## Function : Define MIP fastq file name formats matching regexp
## Returns  : %mip_file_name_regexp
## Arguments:

    my ($arg_href) = @_;

    my %mip_file_name_regexp = (
        features => [qw{ lane date flowcell infile_sample_id index direction }],
        regexp   => q?(\d+)_(\d+)_([^_]+)_([^_]+)_([^_]+)_(\d).fastq?,
    );

    return %mip_file_name_regexp;
}

1;
