package MIP::Cli::Mip::Analyse::Cancer;

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
use MooseX::Types::Moose qw{ Str Int HashRef Num Bool ArrayRef };
use Moose::Util::TypeConstraints;

## MIPs lib
use MIP::Main::Analyse qw{ mip_analyse };

our $VERSION = 0.01;

extends(qw{ MIP::Cli::Mip::Analyse });

command_short_description(q{Cancer analysis});

command_long_description(q{Cancer analysis on panel, wes or wgs sequence data});

command_usage(q{mip <analyse> <cancer> <family_id> --config <config_file>});

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
    use MIP::Get::Analysis qw{ print_program };

    ## Mip analyse cancer parameters
    my @definition_files = (
        catfile( $Bin, qw{ definitions mip_parameters.yaml } ),
        catfile( $Bin, qw{ definitions cancer_parameters.yaml } )
    );

    ## Non mandatory parameter definition keys to check
    my $non_mandatory_parameter_keys_path =
      catfile( $Bin, qw{ definitions non_mandatory_parameter_keys.yaml } );

    ## Mandatory parameter definition keys to check
    my $mandatory_parameter_keys_path =
      catfile( $Bin, qw{ definitions mandatory_parameter_keys.yaml } );

    ### %parameter holds all defined parameters for MIP
    ### mip analyse cancer
    my %parameter;
    foreach my $definition_file (@definition_files) {

        %parameter = (
            %parameter,
            parse_definition_file(
                {
                    define_parameters_path => $definition_file,
                    non_mandatory_parameter_keys_path =>
                      $non_mandatory_parameter_keys_path,
                    mandatory_parameter_keys_path =>
                      $mandatory_parameter_keys_path,
                }
            )
        );
    }

    ## Print programs and exit
    if ( $active_parameter{print_programs} ) {

        print_program(
            {
                define_parameters_files_ref => \@definition_files,
                parameter_href              => \%parameter,
                print_program_mode => $active_parameter{print_program_mode},
            }
        );
        exit;
    }

    ### To add/write parameters in the correct order
    ## Adds the order of first level keys from yaml file to array
    my @order_parameters;
    foreach my $define_parameters_file (@definition_files) {

        push @order_parameters,
          order_parameter_names(
            {
                file_path => $define_parameters_file,
            }
          );
    }

    ## File info hash
    my %file_info = (

        # BWA human genome reference file endings
        bwa_build_reference => [qw{ .bwt .ann .amb .pac .sa }],

        exome_target_bed =>
          [qw{ .infile_list .pad100.infile_list .pad100.interval_list }],

        # Human genome meta files
        human_genome_reference_file_endings => [qw{ .dict .fai }],

    );

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
            cmd_tags    => [q{gatk_baserecalibration_known_sites}],
            documentation =>
              q{Set the references to be decomposed and normalized},
            is  => q{rw},
            isa => ArrayRef [Str],
        )
    );

    option(
        q{exome_target_bed} => (
            cmd_aliases => [qw{ extb }],
            cmd_tags =>
              [q{file.bed=Sample_id; Default: latest_supported_capturekit.bed}],
            documentation => q{Exome target bed file per sample id},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{sample_origin} => (
            cmd_aliases   => [qw{ sao }],
            cmd_tags      => [q{sample_id=sample_origin}],
            documentation => q{Sample origin for analysis},
            is            => q{rw},
            isa           => HashRef,
        )
    );

    option(
        q{pvardict} => (
            cmd_aliases   => [qw{ pvrd }],
            cmd_tags      => [q{Analysis recipe switch}],
            documentation => q{Variant calling using Vardict},
            is            => q{rw},
            isa           => enum( [ 0, 1, 2 ] ),
        )
    );

    option(
        q{vrd_af_threshold} => (
            cmd_aliases   => [qw{ vdraf }],
            cmd_tags      => [q{Default: 0.01}],
            documentation => q{AF threshold for variant calling},
            is            => q{rw},
            isa           => Num,
        )
    );

    option(
        q{vrd_chrom_start} => (
            cmd_aliases   => [qw{ vrdcs }],
            cmd_tags      => [q{Default: 1}],
            documentation => q{Column for chromosome in the output},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{vrd_input_bed_file} => (
            cmd_aliases   => [qw{ vrdbed }],
            documentation => q{Infile path for region info bed file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{vrd_max_mm} => (
            cmd_aliases   => [qw{ vrdmm }],
            cmd_tags      => [q{Default: 4.5}],
            documentation => q{Maximum mean mismatches allowed},
            is            => q{rw},
            isa           => Num,
        )
    );

    option(
        q{vrd_max_pval} => (
            cmd_aliases   => [qw{ vrdmp }],
            cmd_tags      => [q{Default: 0.9}],
            documentation => q{Maximum p-value, set to 0 to keep all variants},
            is            => q{rw},
            isa           => Num,
        )
    );

    option(
        q{vrd_region_end} => (
            cmd_aliases   => [qw{ vrdre }],
            cmd_tags      => [q{Default: 3}],
            documentation => q{Column for region end position in the output},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{vrd_region_start} => (
            cmd_aliases   => [qw{ vrdrs }],
            cmd_tags      => [q{Default: 2}],
            documentation => q{Column for region start position in the output},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{vrd_segment_annotn} => (
            cmd_aliases   => [qw{ vrdsa }],
            cmd_tags      => [q{Default: 4}],
            documentation => q{Column for segment annotation in the output},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{vrd_somatic_only} => (
            cmd_aliases   => [qw{ vrdso }],
            cmd_tags      => [q{Default: no}],
            documentation => q{Output only candidate somatic},
            is            => q{rw},
            isa           => Str,
        )
    );

    return;
}

1;
