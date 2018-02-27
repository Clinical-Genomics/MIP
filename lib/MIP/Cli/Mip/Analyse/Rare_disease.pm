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
use Moose::Util::TypeConstraints;

## MIPs lib
use MIP::Main::Analyse qw{ mip_analyse };

our $VERSION = 0.01;

extends(qw{ MIP::Cli::Mip::Analyse });

command_short_description(q{Rare disease analysis});

command_long_description(
    q{Rare disease analysis on wes, wgs or mixed sequence data});

command_usage(q{mip <analyse> <rare_disease> -pbwa_mem INT});

## Define, check and get Cli supplied parameters
_build_usage();

sub run {
    my ($arg_href) = @_;

    ## Remove Moose::App extra variable
    delete $arg_href->{extra_argv};

    ## Input from Cli
    my %active_parameter = %{$arg_href};

    use MIP::File::Format::Parameter qw{ parse_definition_file  };
    use MIP::File::Format::Yaml qw{ order_parameter_names };

    ## Mip analyse rare_disease parameters
    my $define_parameters_path =
      catfile( $Bin, qw{ definitions rare_disease_parameters.yaml } );

    ## Non mandatory parameter definition keys to check
    my $non_mandatory_parameter_keys_path = catfile( $Bin,
        qw{ definitions rare_disease_non_mandatory_parameter_keys.yaml } );
    ## Mandatory parameter definition keys to check
    my $mandatory_parameter_keys_path = catfile( $Bin,
        qw{ definitions rare_disease_mandatory_parameter_keys.yaml } );

    ### %parameter holds all defined parameters for MIP
    ### analyse rare_disease
    my %parameter = parse_definition_file(
        {
            define_parameters_path => $define_parameters_path,
            non_mandatory_parameter_keys_path =>
              $non_mandatory_parameter_keys_path,
            mandatory_parameter_keys_path => $mandatory_parameter_keys_path,
        }
    );

    ### To add/write parameters in the correct order
    ## Adds the order of first level keys from yaml file to array
    my @order_parameters = order_parameter_names(
        {
            file_path => $define_parameters_path,
        }
    );

    ## File info hash
    my %file_info = (

        # BWA human genome reference file endings
        bwa_build_reference => [qw{ .bwt .ann .amb .pac .sa }],

        exome_target_bed =>
          [qw{ .infile_list .pad100.infile_list .pad100.interval_list }],

        # Human genome meta files
        human_genome_reference_file_endings => [qw{ .dict .fai }],

        # RTG human genome reference file endings
        rtg_vcfeval_reference_genome => [qw{ _sdf_dir }],
    );

    # write_args(\%parameter);
    #write_args( \%active_parameter );

    #exit;

    mip_analyse(
        {
            active_parameter_href => \%active_parameter,
            file_info_href        => \%file_info,
            parameter_href        => \%parameter,
            order_parameters_ref  => \@order_parameters,
        }
    );

    return;
}

sub _build_usage {

## Function : Get and/or set input parameters
## Returns  :
## Arguments:

    option(
        q{decompose_normalize_references} => (
            cmd_aliases => [qw{ dnr }],
            cmd_flag    => q{dec_norm_ref},
            documentation =>
              q{Set the references to be decomposed and normalized},
            is  => q{rw},
            isa => q{ArrayRef},
        )
    );

    option(
        q{expected_coverage} => (
            cmd_aliases   => [qw{ ec }],
            cmd_tags      => [q{sample_id=expected_coverage}],
            documentation => q{Expected mean target coverage for analysis},
            is            => q{rw},
            isa           => q{HashRef},
        )
    );

    option(
        q{pbwa_mem} => (
            cmd_aliases   => [qw{ pmem }],
            cmd_flag      => q{pbwa_mem},
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Align reads using Bwa Mem},
            is            => q{rw},
            isa           => enum( [ 1, 2 ] ),
        )
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
