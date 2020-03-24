#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2018 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.03;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::File_info}      => [qw{ check_parameter_metafiles }],
        q{MIP::Io::Read}       => [qw{ read_from_file }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File_info qw{ check_parameter_metafiles };
use MIP::Io::Read qw{ read_from_file };

diag(   q{Test check_parameter_metafiles from File_info.pm v}
      . $MIP::File_info::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my %parameter = read_from_file(
    {
        format => q{yaml},
        path   => catfile( dirname($Bin), qw{ definitions rd_dna_parameters.yaml} ),
    }
);

## File info hash
my %file_info = (

    # BWA human genome reference file endings
    bwa_build_reference => [qw{ .bwt .ann .amb .pac .sa }],

    exome_target_bed => [qw{ .interval_list .pad100.interval_list }],

    # Human genome meta files
    human_genome_reference_file_endings => [qw{ .dict .fai }],

    # RTG human genome reference file endings
    rtg_vcfeval_reference_genome => [qw{ _sdf_dir }],
);

my $parameter_name = q{exome_target_bed};

## Given hash entries with no active parameter and existing files
my %active_parameter = (
    not_correct_key => {
        catfile( $Bin,
            qw{ data references grch37_agilent_sureselect_targets_cre_-v1-.bed } ) =>
          q{sample1},
    },
);

check_parameter_metafiles(
    {
        active_parameter_href => \%active_parameter,
        file_info_href        => \%file_info,
        parameter_href        => \%parameter,
    }
);

## Then set build switch should not be set to zero i.e. default of "1" is kept
## as this will not be evaluated in sub or downstream
is( $parameter{$parameter_name}{build_file}, 1,
    q{No active parameter and existing file} );

## Given hash entries with active parameter, existing files, and no active associated programs
%active_parameter = (
    exome_target_bed => {
        catfile( $Bin,
            qw{ data references grch37_agilent_sureselect_targets_cre_-v1-.bed } ) =>
          q{sample1},
    },
);

check_parameter_metafiles(
    {
        active_parameter_href => \%active_parameter,
        file_info_href        => \%file_info,
        parameter_href        => \%parameter,
    }
);

## Then set build switch should not be set to zero i.e. default of "1" is kept
## as this will not be evaluated in sub or downstream
is( $parameter{$parameter_name}{build_file},
    1, q{Active parameter and existing file, and no active associated program} );

## Given hash entries with active parameter, files exists, and active associated programs
%active_parameter = (
    exome_target_bed => {
        catfile( $Bin,
            qw{ data references grch37_agilent_sureselect_targets_cre_-v1-.bed } ) =>
          q{sample1},
    },
    picardtools_collecthsmetrics => 1,
);

check_parameter_metafiles(
    {
        active_parameter_href => \%active_parameter,
        file_info_href        => \%file_info,
        parameter_href        => \%parameter,
    }
);

## Then set build switch should be set to zero to not build meta files as
## they are all ready and do not need to be built
is( $parameter{$parameter_name}{build_file},
    0, q{Active parameter, existing file and active associated program} );

## Given hash entries where one path exists and one does not
$active_parameter{exome_target_bed}
  { catfile( $Bin, qw{ data references does_not_exists.bed } ) } = q{sample1};

check_parameter_metafiles(
    {
        active_parameter_href => \%active_parameter,
        file_info_href        => \%file_info,
        parameter_href        => \%parameter,
    }
);

## Then set build file to true to rebuild for all parameter keys
is( $parameter{exome_target_bed}{build_file},
    1, q{Set build file switch for hash parameter reference with mixed existence to 1} );

## Given scalar entries with active parameter, files exists, and active associated programs
$active_parameter{bwa_build_reference} = 1;
$active_parameter{human_genome_reference} =
  catfile( $Bin, qw{ data references grch37_homo_sapiens_-d5-.fasta } );
$active_parameter{bwa_mem} = 1;

check_parameter_metafiles(
    {
        active_parameter_href => \%active_parameter,
        file_info_href        => \%file_info,
        parameter_href        => \%parameter,
    }
);

## Then set build file to false
is( $parameter{bwa_build_reference}{build_file}, 0,
q{Active parameter, existing file and active associated program for scalar parameter reference}
);

## Given scalar entries with active parameter, files do not exists, and active associated programs
%file_info = (

    # BWA human genome reference file endings
    bwa_build_reference => [qw{ .does_not_exist .ann .amb .pac .sa }],
);

check_parameter_metafiles(
    {
        active_parameter_href => \%active_parameter,
        file_info_href        => \%file_info,
        parameter_href        => \%parameter,
    }
);

## Then set build file to true
is( $parameter{bwa_build_reference}{build_file}, 1,
q{Active parameter, no existing file and active associated program for scalar parameter reference}
);

done_testing();
