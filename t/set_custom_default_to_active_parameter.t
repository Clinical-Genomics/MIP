#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
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
use MIP::Constants qw { $COMMA $SPACE $UNDERSCORE };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.11;

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
        q{MIP::File::Format::Yaml} => [qw{ load_yaml }],
        q{MIP::Get::Parameter}     => [qw{ get_capture_kit }],
        q{MIP::Set::Parameter}     => [qw{ set_custom_default_to_active_parameter }],
        q{MIP::Test::Fixtures}     => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File::Format::Yaml qw{ load_yaml };
use MIP::Get::Parameter qw{ get_capture_kit };
use MIP::Set::Parameter qw{ set_custom_default_to_active_parameter };

diag(   q{Test set_custom_default_to_active_parameter from Parameter.pm v}
      . $MIP::Set::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Creates log object
my $log = test_log( { no_screen => 1 } );

my %active_parameter = (
    cluster_constant_path  => catfile(qw{ constant path }),
    conda_path             => catdir( $Bin, qw{ data modules miniconda } ),
    case_id                => 1,
    human_genome_reference => catfile(qw{ a test grch37_human_genom_reference.fasta }),
    load_env               => {
        test_env => {
            gatk                   => undef,
            method                 => q{conda},
            mip                    => undef,
            varianteffectpredictor => undef,
        },
        test_env_1 => {
            method => q{conda},
            picard => undef,
        },
    },
    outdata_dir     => catfile(qw{ a outdata dir }),
    sample_ids      => [qw{ sample_1 }],
    select_programs => undef,
);

## Mip analyse rd_dna parameters
## The order of files in @definition_files should follow commands inheritance
my @definition_files = (
    catfile( dirname($Bin), qw{ definitions analyse_parameters.yaml } ),
    catfile( dirname($Bin), qw{ definitions rd_dna_parameters.yaml } ),
);

my %parameter;

DEFINITION_FILE:
foreach my $definition_file (@definition_files) {

    %parameter = (
        %parameter,
        load_yaml(
            {
                yaml_file => $definition_file,
            }
        ),
    );
}

## Given custom default parameters
my @custom_default_parameters = qw{ analysis_type
  bwa_build_reference
  exome_target_bed
  infile_dirs
  reference_dir
  rtg_vcfeval_reference_genome
  sample_info_file
  select_programs
  sv_vcfparser_select_file
  temp_directory
  vcfparser_select_file };

PARAMETER_NAME:
foreach my $parameter_name (@custom_default_parameters) {

    set_custom_default_to_active_parameter(
        {
            active_parameter_href => \%active_parameter,
            parameter_href        => \%parameter,
            parameter_name        => $parameter_name,
        }
    );
}

my $vcfparser_select_file_path = catfile(
    $active_parameter{cluster_constant_path},
    $active_parameter{case_id},
    q{gene_panel.bed}
);

my %expected_default = (
    bwa_build_reference => {
        default    => $active_parameter{human_genome_reference},
        test_label => q{Set human_genome_reference default for bwa},
    },
    reference_dir => {
        default    => cwd(),
        test_label => q{Set reference_dir default },
    },
    rtg_vcfeval_reference_genome => {
        default => $active_parameter{human_genome_reference},
        test_label =>
          q{Set human_genome_reference default for rtg vcfeval reference genome},
    },
    sv_vcfparser_select_file => {
        default    => $vcfparser_select_file_path,
        test_label => q{Set sv_vcfparser_select_file default },
    },
    temp_directory => {
        default    => catfile( $active_parameter{outdata_dir}, q{$SLURM_JOB_ID} ),
        test_label => q{Set temp_directory default },
    },
    vcfparser_select_file => {
        default    => $vcfparser_select_file_path,
        test_label => q{Set vcfparser_select_file default },
    },
);

## Then the defaults should be set for each parameter given
is( $active_parameter{analysis_type}{sample_1}, q{wgs}, q{Set analysis_type default} );

while ( my ( $parameter_name, $meta_data_href ) = each %expected_default ) {

    is(
        $active_parameter{$parameter_name},
        $meta_data_href->{default},
        $meta_data_href->{test_label}
    );
}

## Given an download pipe
$active_parameter{download_pipeline_type} = q{rd_dna};

set_custom_default_to_active_parameter(
    {
        active_parameter_href => \%active_parameter,
        parameter_href        => \%parameter,
        parameter_name        => q{temp_directory},
    }
);

is(
    $active_parameter{temp_directory},
    catfile( cwd(), qw{mip_download $SLURM_JOB_ID} ),
    q{Set temp_directory default for download pipe }
);

## Given an analysis type, when unset for a sample
# Clear analysis type
delete $active_parameter{analysis_type}{sample_1};

set_custom_default_to_active_parameter(
    {
        active_parameter_href => \%active_parameter,
        parameter_href        => \%parameter,
        parameter_name        => q{infile_dirs},
    }
);

## Then the defaults should be set for analysis type before setting default for infile dirs
is( $active_parameter{analysis_type}{sample_1},
    q{wgs}, q{Set analysis_type default within infile dirs} );

CAPTURE_KIT:
foreach my $capture_kit ( keys %{ $active_parameter{exome_target_bed} } ) {

    is(
        $capture_kit,
        $parameter{supported_capture_kit}{default}{latest},
        q{Set default capture kit}
    );

    is( $active_parameter{exome_target_bed}{$capture_kit},
        q{sample_1}, q{Set default capture kit for sample_1} );

}

my $sample_info_file = catfile(
    $active_parameter{outdata_dir},
    $active_parameter{case_id},
    $active_parameter{case_id} . $UNDERSCORE . q{qc_sample_info.yaml}
);

is( $parameter{sample_info_file}{default},
    $sample_info_file, q{Set default sample_info_file} );

is( $parameter{qccollect_sampleinfo_file}{default},
    $sample_info_file, q{Set default qccollect sample_info_file} );

my $path = catfile(
    $active_parameter{cluster_constant_path},
    $active_parameter{case_id},
    $active_parameter{analysis_type}{sample_1},
    q{sample_1}, q{fastq}
);
is( $active_parameter{infile_dirs}{$path}, q{sample_1}, q{Set default infile_dirs} );

## Test setting custom paths
my %test_hash = (
    gatk_path        => catdir( $Bin, qw{ data modules GenomeAnalysisTK-3.7 } ),
    picardtools_path => catdir(
        $active_parameter{conda_path}, qw{ envs test_env_1 share picard-2.14.1-0 }
    ),
    snpeff_path =>
      catdir( $active_parameter{conda_path}, qw{ envs test_env share snpeff } ),
    vep_directory_path =>
      catdir( $active_parameter{conda_path}, qw{ envs test_env ensembl-vep } ),
);

TEST_PATH:
foreach my $test_path ( keys %test_hash ) {

    set_custom_default_to_active_parameter(
        {
            active_parameter_href => \%active_parameter,
            parameter_href        => \%parameter,
            parameter_name        => $test_path,
        }
    );

    is( $active_parameter{$test_path},
        $test_hash{$test_path}, q{Set default} . $SPACE . $test_path );
}

done_testing();
