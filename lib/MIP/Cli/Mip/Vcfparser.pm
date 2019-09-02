package MIP::Cli::Mip::Vcfparser;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

$Params::Check::PRESERVE_CASE = 1;    #Do not convert to lower case

## CPANM
use autodie qw{ :all };
use Modern::Perl qw{ 2018 };
use Moose::Util::TypeConstraints;
use MooseX::App::Command;
use MooseX::Types::Moose qw{ ArrayRef Bool HashRef Int Str };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ %ANALYSIS $COLON $DASH $NEWLINE };
use MIP::Log::MIP_log4perl qw{ initiate_logger };
use MIP::Main::Vcfparser qw{ mip_vcfparser };

our $VERSION = 1.00;

command_short_description(q{MIP vcfparser command});

command_long_description(
    q{Entry point for splitting VCF into clinical and research variants});

command_usage(q{vcfparser [options] infile.vcf [OPTIONS] > outfile.vcf});

## Constants
Readonly my $ANNOTATION_DISTANCE => $ANALYSIS{ANNOTATION_DISTANCE};

## Define, check and get Cli supplied parameters
_build_usage();

sub run {

    ## Input from Cli
    my ($arg_href) = @_;

    ## Remove Moose::App extra variable
    delete $arg_href->{extra_argv};

    # Flatten argument(s)
    my $parse_vep          = $arg_href->{parse_vep};
    my $range_feature_file = $arg_href->{range_feature_file};
    my @range_feature_annotation_columns;
    if ( exists $arg_href->{range_feature_annotation_columns} ) {

        @range_feature_annotation_columns =
          @{ $arg_href->{range_feature_annotation_columns} };
    }
    my $select_feature_file            = $arg_href->{select_feature_file};
    my $select_feature_matching_column = $arg_href->{select_feature_matching_column};
    my @select_feature_annotation_columns;
    if ( exists $arg_href->{select_feature_annotation_columns} ) {

        @select_feature_annotation_columns =
          @{ $arg_href->{select_feature_annotation_columns} };
    }
    my $select_outfile       = $arg_href->{select_outfile};
    my $padding              = $arg_href->{padding};
    my $per_gene             = $arg_href->{per_gene};
    my $pli_values_file_path = $arg_href->{pli_values_file};
    my $write_software_tag   = $arg_href->{write_software_tag};
    my $log_file             = $arg_href->{log_file};

    ## STDIN
    my $infile = $arg_href->{infile};

    # Create anonymous filehandle
    my $VCF_IN_FH = IO::Handle->new();

    ## Enables cmd "vcfparser" to print usage help
    if ( $infile eq $DASH ) {

        $VCF_IN_FH = *STDIN;
    }
    else {
        open $VCF_IN_FH, q{<}, $infile
          or croak( q{Cannot open } . $infile . $COLON . $OS_ERROR, $NEWLINE );
    }

    ## Creates log object
    my $log = initiate_logger(
        {
            file_path => $log_file,
            log_name  => q{Vcfparser},
        }
    );

    ## Basic flag option check
    if ( not @range_feature_annotation_columns and $range_feature_file ) {

        $log->fatal(
            q{Need to specify which feature column(s) to use with range feature file: }
              . $range_feature_file
              . q{ when annotating variants by using flag -rf_ac},
            $NEWLINE
        );
        exit 1;
    }
    if ( not $select_feature_matching_column and $select_feature_file ) {

        $log->fatal(
            q{Need to specify which feature column to use with select feature file: }
              . $select_feature_file
              . q{ when selecting variants by using flag -sf_mc},
            $NEWLINE
        );
        exit 1;
    }
    if ( not $select_outfile and $select_feature_file ) {

        $log->fatal(
q{Need to specify which a select outfile to use when selecting variants by using flag -sof},
            $NEWLINE
        );
        exit 1;
    }

    mip_vcfparser(
        {
            VCF_IN_FH                            => $VCF_IN_FH,
            padding                              => $padding,
            parse_vep                            => $parse_vep,
            per_gene                             => $per_gene,
            pli_values_file_path                 => $pli_values_file_path,
            range_feature_annotation_columns_ref => \@range_feature_annotation_columns,
            range_feature_file                   => $range_feature_file,
            select_feature_file                  => $select_feature_file,
            select_feature_matching_column       => $select_feature_matching_column,
            select_outfile_path                  => $select_outfile,
            write_software_tag                   => $write_software_tag,
        }
    );

    return;
}

sub _build_usage {

## Function : Get and/or set input parameters
## Returns  :
## Arguments:

    my ($arg_href) = @_;

    parameter(
        q{infile} => (
            documentation => q{Infile path},
            is            => q{rw},
            isa           => Str,
            required      => 1,
        )
    );

    option(
        q{parse_vep} => (
            cmd_aliases   => [qw{ pvep }],
            documentation => q{Parse VEP transcript specific entries},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{log_file} => (
            cmd_aliases   => [qw{ l }],
            default       => catfile( cwd(), q{vcfparser.log} ),
            documentation => q{Log file},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{padding} => (
            cmd_aliases   => [qw{ pad }],
            cmd_flag      => q{padding},
            default       => $ANNOTATION_DISTANCE,
            documentation => q{Number of nucleotides to pad},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{per_gene} => (
            cmd_aliases   => [qw{ peg }],
            documentation => q{Output most severe annotations},
            is            => q{rw},
            isa           => Bool,
        )
    );

    option(
        q{pli_values_file} => (
            cmd_aliases   => [qw{ pli }],
            cmd_tags      => [q{Format: TSV}],
            documentation => q{Pli value file path},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{range_feature_file} => (
            cmd_aliases   => [qw{ rf }],
            cmd_tags      => [q{Format: TSV}],
            documentation => q{Range feature file path},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{range_feature_annotation_columns} => (
            cmd_aliases   => [qw{ rf_ac }],
            cmd_flag      => q{range_feature_annotation_columns},
            documentation => q{Range feature file annotation columns},
            is            => q{rw},
            isa           => ArrayRef,
        )
    );

    option(
        q{select_feature_file} => (
            cmd_aliases   => [qw{ sf }],
            cmd_tags      => [q{Format: TSV}],
            documentation => q{Select feature file path},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{select_feature_matching_column} => (
            cmd_aliases   => [qw{ sf_mc }],
            cmd_flag      => q{select_feature_matching_column},
            documentation => q{Select feature file annotation columns},
            is            => q{rw},
            isa           => Int,
        )
    );

    option(
        q{select_feature_annotation_columns} => (
            cmd_aliases   => [qw{ sf_ac }],
            cmd_flag      => q{select_feature_annotation_columns},
            documentation => q{Select feature file annotation columns},
            is            => q{rw},
            isa           => ArrayRef,
        )
    );

    option(
        q{select_outfile} => (
            cmd_aliases   => [qw{ sof }],
            cmd_tags      => [q{Format: VCF}],
            documentation => q{Select feature outfile},
            is            => q{rw},
            isa           => Str,
        )
    );

    option(
        q{write_software_tag} => (
            cmd_aliases   => [qw{ wst }],
            default       => 1,
            documentation => q{Write software tag},
            is            => q{rw},
            isa           => Bool,
        )
    );

    return;
}

1;
