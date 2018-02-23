package MIP::Cli::Mip::Analyse::Rare_disease;

use Carp;
use File::Spec::Functions qw{ catfile };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use MooseX::App::Command;

## MIPs lib
use MIP::Check::Parameter qw{  check_parameter_hash };
use MIP::File::Format::Yaml qw{ load_yaml order_parameter_names };
use MIP::Scripts::Analyse qw{ mip_analyse };

our $VERSION = 0.01;

extends(qw{ MIP::Cli::Mip::Analyse });

command_short_description(q{Rare disease analysis});

command_long_description(
    q{Rare disease analysis on wes, wgs or mixed sequence data});

command_usage(q{mip <analyse> <rare_disease> -pbwa_mem INT});

option(
    q{pbwa_mem} => (
        cmd_aliases   => [qw{ pmem }],
        cmd_flag      => q{pbwa_mem},
        cmd_tags      => [q{Analysis recipe switch}],
        documentation => q{Align reads using Bwa Mem (defaults to "0" (=no))},
        is            => q{rw},
        isa           => q{Int},
        required      => 0,
    )
);

sub run {
    my ($arg_href) = @_;

    ## Remove Moose::App extra variable
    delete $arg_href->{extra_argv};

    ## Input from Cli
    my %active_parameter = %{$arg_href};

    ## Load definition file of mip analyse rare_disease parameters
    my $definitions_file =
      catfile( $Bin, qw{ definitions rare_disease_parameters.yaml } );

    ### %parameter holds all defined parameters for MIP
    ### analyse rare_disease

    #########
    ### TO DO: Add parse_definition_file sub for order and check return aprameter
    #########

    ## Loads a YAML file into an arbitrary hash and returns it.
    my %parameter = load_yaml( { yaml_file => $definitions_file, } );

    ### To add/write parameters in the correct order
    ## Adds the order of first level keys from yaml file to array
    my @order_parameters = order_parameter_names(
        {
            file_path => $definitions_file,
        }
    );

## Load mandatory keys and values for parameters
    my %mandatory_key = load_yaml(
        {
            yaml_file => catfile(
                $Bin,
                qw{ definitions rare_disease_mandatory_parameter_keys.yaml }
            ),
        }
    );

## Load non mandatory keys and values for parameters
    my %non_mandatory_key = load_yaml(
        {
            yaml_file => catfile(
                $Bin,
                qw{ definitions rare_disease_non_mandatory_parameter_keys.yaml }
            ),

        }
    );

    ## Eval parameter hash
    check_parameter_hash(
        {
            file_path              => $definitions_file,
            non_mandatory_key_href => \%non_mandatory_key,
            mandatory_key_href     => \%mandatory_key,
            parameter_href         => \%parameter,
        }
    );

    ### To DO new sub ends here with return of parameter hash

    # write_args(\%parameter);
    write_args( \%active_parameter );

    #exit;

    mip_analyse(
        {
            active_parameter_href => \%active_parameter,
            parameter_href        => \%parameter,
            order_parameters_ref  => \@order_parameters,
        }
    );

    return;
}

sub write_args {

    my ($arg_href) = @_;

    # do something
    use Data::Dumper;

    #    say STDERR $arg_href->{pbwa_mem};
    #    foreach my $sample ( @{ $arg_href->{sample_ids} } ) {
    #        say STDERR $sample;
    #    }
    print Dumper($arg_href);
    return;
}

1;
