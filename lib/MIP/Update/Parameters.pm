package MIP::Update::Parameters;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use Readonly;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ update_dynamic_config_parameters update_reference_parameters update_vcfparser_outfile_counter };
}

## Constants
Readonly my $SPACE => q{ };

sub update_dynamic_config_parameters {

## Function : Updates the config file to particular user/cluster for dynamic config parameters following specifications. Leaves other entries untouched.
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $parameter_name        => MIP Parameter to update

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $parameter_name;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        parameter_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$parameter_name
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return if ( not defined $active_parameter_href->{$parameter_name} );

    my @dynamic_parameters =
      qw{ cluster_constant_path analysis_constant_path family_id outaligner_dir };

  DYNAMIC_PARAMETER:
    foreach my $dynamic_parameter (@dynamic_parameters) {

        ## Replace dynamic config parameters with actual value that is now set from cmd or config
        $active_parameter_href->{$parameter_name} =~
s/$dynamic_parameter!/$active_parameter_href->{$dynamic_parameter}/smgi;
    }
    return;
}

sub update_reference_parameters {

## Function : Update reference parameters with mip_reference directory path
## Returns  :
## Arguments: $active_parameter_href   => Holds all set parameter for analysis
##          : $associated_programs_ref => The parameters program(s) {REF}
##          : $parameter_name          => Parameter name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $associated_programs_ref;
    my $parameter_name;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
        associated_programs_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$associated_programs_ref
        },
        parameter_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$parameter_name,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Set::Parameter qw{ set_parameter_reference_dir_path };

    ## Check all programs that use parameter
  ASSOCIATED_PROGRAM:
    foreach my $associated_program ( @{$associated_programs_ref} ) {

        my $program_name = $active_parameter_href->{$associated_program};

        ## Only check active programs parameters
        next ASSOCIATED_PROGRAM if ( not $program_name );

        ## Update path for supplied reference(s) associated with
        ## parameter that should reside in the mip reference directory to full path
        set_parameter_reference_dir_path(
            {
                active_parameter_href => $active_parameter_href,
                parameter_name        => $parameter_name,
            }
        );
        ## Only need to perform update once per parameter
        return;
    }
    return;
}

sub update_vcfparser_outfile_counter {

## Function : Determine the number of outfile after vcfparser
## Returns  :
## Arguments: $active_parameter_href   => Holds all set parameter for analysis

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;

    my $tmpl = {
        active_parameter_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$active_parameter_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Create link
    my %vcfparser_select_file = (
        pvcfparser => { vcfparser_select_file => q{vcfparser_outfile_count} },
        psv_vcfparser =>
          { sv_vcfparser_select_file => q{sv_vcfparser_outfile_count} },
    );
## Determine if to expect select outfile for vcfparser and sv_vcfparser
  PROGRAM:
    foreach my $program ( keys %vcfparser_select_file ) {

        next PROGRAM if ( not $active_parameter_href->{$program} );

      FILES:
        while ( my ( $parameter_name, $parameter_name_counter ) =
            each %{ $vcfparser_select_file{$program} } )
        {

            $active_parameter_href->{$parameter_name_counter} =
              _set_vcfparser_file_counter(
                {
                    parameter_name => $active_parameter_href->{$parameter_name},
                }
              );
        }
    }
    return;
}

sub _set_vcfparser_file_counter {

## Function : Set the expected number of outputfile after vcfparser
## Returns  :
## Arguments: $parameter_name => Vcfparser select file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_name;

    my $tmpl =
      { parameter_name =>
          { required => 1, strict_type => 1, store => \$parameter_name }, };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

## To track if vcfparser was used with a vcfparser_select_file (=2) or not (=1)
    if ( not defined $parameter_name ) {

        ## No select file was given
        return 1;
    }
    else {

        ## To track if vcfparser was used with a vcfparser_select_file (=2) or not (=1)
        return 2;
    }
    return;
}

1;
