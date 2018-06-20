package MIP::Get::Analysis;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## Third party module(s)
use autodie;
use List::MoreUtils qw{ all any };
use Log::Log4perl;
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.05;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ get_dependency_tree get_dependency_tree_chain get_dependency_tree_order get_overall_analysis_type print_program };

}

## Constants
Readonly my $SPACE   => q{ };
Readonly my $NEWLINE => qq{\n};

sub get_dependency_tree {

## Function  : Collects all downstream programs from initation point.
## Returns   :
## Arguments : $current_chain           => Current chain
##           : $is_program_found_ref    => Found initiation program {REF}
##           : $is_chain_found_ref      => Found program chain
##           : $program                 => Initiation point
##           : $start_with_programs_ref => Store programs
##           : $dependency_tree_href    => Dependency hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $current_chain;
    my $is_program_found_ref;
    my $is_chain_found_ref;
    my $program;
    my $start_with_programs_ref;
    my $dependency_tree_href;

    my $tmpl = {
        current_chain => {
            store       => \$current_chain,
            strict_type => 1,
        },
        is_program_found_ref => {
            default     => \$$,
            store       => \$is_program_found_ref,
            strict_type => 1,
        },
        is_chain_found_ref => {
            default     => \$$,
            store       => \$is_chain_found_ref,
            strict_type => 1,
        },
        program => {
            required    => 1,
            store       => \$program,
            strict_type => 1,
        },
        start_with_programs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$start_with_programs_ref,
            strict_type => 1,
        },
        dependency_tree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$dependency_tree_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Copy hash to enable recursive removal of keys
    my %tree = %{$dependency_tree_href};

  KEY_VALUE_PAIR:
    while ( my ( $key, $value ) = each %tree ) {

        ## Do not enter into more chains than one if program and chain is found
        next KEY_VALUE_PAIR
          if ( $key =~ /CHAIN_/sxm
            && ${$is_program_found_ref}
            && ${$is_chain_found_ref} ne q{CHAIN_MAIN} );

        ## Do not add program name or PARALLEL
        if ( $key =~ /CHAIN_/sxm ) {

            $current_chain = $key;
        }
        ## Call recursive
        if ( ref $value eq q{HASH} ) {

            get_dependency_tree(
                {
                    current_chain           => $current_chain,
                    dependency_tree_href    => $value,
                    is_program_found_ref    => $is_program_found_ref,
                    is_chain_found_ref      => $is_chain_found_ref,
                    program                 => $program,
                    start_with_programs_ref => $start_with_programs_ref,
                }
            );
        }
        elsif ( ref $value eq q{ARRAY} ) {
            ## Inspect element

          ELEMENT:
            foreach my $element ( @{$value} ) {

                ## Call recursive
                if ( ref $element eq q{HASH} ) {

                    get_dependency_tree(
                        {
                            current_chain           => $current_chain,
                            dependency_tree_href    => $element,
                            is_program_found_ref    => $is_program_found_ref,
                            is_chain_found_ref      => $is_chain_found_ref,
                            program                 => $program,
                            start_with_programs_ref => $start_with_programs_ref,
                        }
                    );
                }
                ## Found initiator program
                if ( ref $element ne q{HASH}
                    && $element eq $program )
                {

                    ## Start collecting programs downstream
                    ${$is_program_found_ref} = 1;

                    ## Found chain that program belongs to
                    # Set is part of chain signal
                    ${$is_chain_found_ref} = $current_chain;

                }

                ## Special case for parallel section
                if ( $key eq q{PARALLEL}
                    && ${$is_program_found_ref} )
                {

                    if ( any { $_ eq $program } @{$value} ) {

                        ## Add only start_with program from parallel section
                        push @{$start_with_programs_ref}, $program;

                        ## Skip any remaining hash_ref or element
                        last ELEMENT;
                    }
                }

                ## Add downstream programs
                if ( ref $element ne q{HASH}
                    && ${$is_program_found_ref} )
                {

                    push @{$start_with_programs_ref}, $element;
                }
            }
        }

        ## Remove identifier
        delete $tree{$key};
    }
    return;
}

sub get_dependency_tree_chain {

## Function  : Sets chain id to parameters hash from the dependency tree
## Returns   :
## Arguments : $current_chain        => Current chain
##           : $dependency_tree_href => Dependency hash {REF}
##           : $parameter_href       => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $current_chain;
    my $dependency_tree_href;
    my $parameter_href;

    my $tmpl = {
        current_chain => {
            store       => \$current_chain,
            strict_type => 1,
        },
        dependency_tree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$dependency_tree_href,
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

    ## Copy hash to enable recursive removal of keys
    my %tree = %{$dependency_tree_href};

  KEY_VALUE_PAIR:
    while ( my ( $key, $value ) = each %tree ) {

        ## Add ID of chain
        my ($chain_id) = $key =~ /CHAIN_(\S+)/sxm;

        ## If chain_id is found
        if ( defined $chain_id ) {

            ## Set current chain
            $current_chain = $chain_id;
        }

        ## Call recursive
        if ( ref $value eq q{HASH} ) {

            get_dependency_tree_chain(
                {
                    current_chain        => $current_chain,
                    dependency_tree_href => $value,
                    parameter_href       => $parameter_href,
                }
            );
        }
        elsif ( ref $value eq q{ARRAY} ) {
            ## Inspect element

          ELEMENT:
            foreach my $element ( @{$value} ) {

                ## Call recursive
                if ( ref $element eq q{HASH} ) {

                    get_dependency_tree_chain(
                        {
                            current_chain        => $current_chain,
                            dependency_tree_href => $element,
                            parameter_href       => $parameter_href,
                        }
                    );
                }

                ## Found programs
                if ( ref $element ne q{HASH} ) {

                    $parameter_href->{$element}{chain} = $current_chain;
                }

                if ( $key eq q{PARALLEL} ) {

                    $parameter_href->{$element}{chain} = uc $element;
                }

                ## Hash in PARALLEL section create anonymous chain ID
                ## E.g. haplotypecaller->genotypegvcfs
                if ( $key eq uc $element ) {

                  PROGRAM:
                    foreach my $program ( @{$value} ) {

                        $parameter_href->{$program}{chain} = uc $element;
                    }
                    last ELEMENT;
                }
            }
        }

        ## Remove identifier
        delete $tree{$key};
    }
    return;
}

sub get_dependency_tree_order {

## Function  : Collects order of all programs from initiation.
## Returns   :
## Arguments : $programs_ref         => Programs {REF}
##           : $dependency_tree_href => Dependency hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $programs_ref;
    my $dependency_tree_href;

    my $tmpl = {
        programs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$programs_ref,
            strict_type => 1,
        },
        dependency_tree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$dependency_tree_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Copy hash to enable recursive removal of keys
    my %tree = %{$dependency_tree_href};

  KEY_VALUE_PAIR:
    while ( my ( $key, $value ) = each %tree ) {

        ## Call recursive
        if ( ref $value eq q{HASH} ) {

            get_dependency_tree_order(
                {
                    dependency_tree_href => $value,
                    programs_ref         => $programs_ref,
                }
            );
        }
        elsif ( ref $value eq q{ARRAY} ) {
            ## Inspect element

          ELEMENT:
            foreach my $element ( @{$value} ) {

                ## Call recursive
                if ( ref $element eq q{HASH} ) {

                    get_dependency_tree_order(
                        {
                            dependency_tree_href => $element,
                            programs_ref         => $programs_ref,
                        }
                    );
                }
                ## Found program
                if ( ref $element ne q{HASH} ) {

                    ## Add to order
                    push @{$programs_ref}, $element;
                }
            }
        }

        ## Remove identifier
        delete $tree{$key};
    }
    return;
}

sub get_overall_analysis_type {

## Function : Detect if all samples has the same sequencing type and return consensus or mixed
## Returns  : q{consensus} | q{mixed} - analysis_type
## Arguments: $analysis_type_href => Analysis_type hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $analysis_type_href;

    my $tmpl = {
        analysis_type_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$analysis_type_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger(q{MIP});

    my @analysis_types = (qw{ wes wgs wts cancer });

  ANALYSIS:
    foreach my $analysis_type (@analysis_types) {

        ## If consensus is reached
        if ( all { $_ eq $analysis_type } values %{$analysis_type_href} ) {

            return $analysis_type;
        }
    }

    ## Check that the user supplied analysis type is supported
    foreach my $user_analysis_type ( values %{$analysis_type_href} ) {

        if ( not any { $_ eq $user_analysis_type } @analysis_types ) {

            $log->fatal( q{'}
                  . $user_analysis_type
                  . q{' is not a supported analysis_type} );
            $log->fatal( q{Supported analysis types are '}
                  . join( q{', '}, @analysis_types )
                  . q(') );
            $log->fatal(q{Aborting run});
            exit 1;
        }
    }

    # No consensus, then it must be mixed
    return q{mixed};
}

sub print_program {

## Function : Print all supported programs in '-ppm' mode
## Returns  :
## Arguments: $define_parameters_files_ref => MIPs define parameters file
##          : $parameter_href         => Parameter hash {REF}
##          : $print_program_mode     => Mode to run modules in

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_href;

    ## Default(s)
    my $define_parameters_files_ref;
    my $print_program_mode;

    my $tmpl = {
        define_parameters_files_ref => {
            default     => [],
            store       => \$define_parameters_files_ref,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        print_program_mode => {
            allow => [ undef, 0, 1, 2 ],
            default => $arg_href->{print_program_mode} //= 2,
            store => \$print_program_mode,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Format::Yaml qw{ order_parameter_names };
    use MIP::Set::Parameter qw{ set_dynamic_parameter };

    my @printed_programs;

    set_dynamic_parameter(
        {
            aggregates_ref => [q{type:program}],
            parameter_href => $parameter_href,
        }
    );

    ## Adds the order of first level keys from yaml file to array
    my @order_parameters;
    foreach my $define_parameters_file ( @{$define_parameters_files_ref} ) {

        push @order_parameters,
          order_parameter_names(
            {
                file_path => $define_parameters_file,
            }
          );
    }

  PARAMETER:
    foreach my $parameter (@order_parameters) {

        ## Only process programs
        if (
            any { $_ eq $parameter }
            @{ $parameter_href->{dynamic_parameter}{program} }
          )
        {

            if (
                not $parameter =~
                / pbamcalibrationblock | pvariantannotationblock /xsm )
            {

                print {*STDOUT} q{--}
                  . $parameter
                  . $SPACE
                  . $print_program_mode
                  . $SPACE;

                push @printed_programs,
                  $parameter . $SPACE . $print_program_mode;
            }
        }
    }
    print {*STDOUT} $NEWLINE;

    return @printed_programs;
}

1;
