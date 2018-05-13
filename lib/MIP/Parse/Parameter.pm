package MIP::Parse::Parameter;

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

BEGIN {

    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ parse_start_with_program };
}

sub parse_start_with_program {

## Function : Get initiation program, downstream dependencies and update program modes fo start_with_program parameter
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $initiation_file       => Initiation file for pipeline
##          : $parameter_href        => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $initiation_file;
    my $parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        initiation_file => {
            defined     => 1,
            required    => 1,
            store       => \$initiation_file,
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

    use MIP::File::Format::Yaml qw{ load_yaml };
    use MIP::Get::Analysis qw{ get_dependency_tree };
    use MIP::Update::Programs qw{  update_program_mode_with_start_with };

    return if ( not defined $active_parameter_href->{start_with_program} );

    my %dependency_tree = load_yaml(
        {
            yaml_file => $initiation_file,
        }
    );

    my @start_with_programs;
    my $is_program_found = 0;
    my $is_chain_found   = 0;

    ## Collects all downstream programs from initation point
    get_dependency_tree(
        {
            dependency_tree_href => \%dependency_tree,
            is_program_found_ref => \$is_program_found,
            is_chain_found_ref   => \$is_chain_found,
            program => $active_parameter_href->{start_with_program},
            start_with_programs_ref => \@start_with_programs,
        }
    );

    ## Update program mode depending on start with flag
    update_program_mode_with_start_with(
        {
            active_parameter_href => $active_parameter_href,
            programs_ref => \@{ $parameter_href->{dynamic_parameter}{program} },
            start_with_programs_ref => \@start_with_programs,
        }
    );
    return 1;
}

1;
