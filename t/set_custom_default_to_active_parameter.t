#!/usr/bin/env perl

use 5.018;
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
use Modern::Perl qw{ 2014 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.2;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $COMMA      => q{,};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

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
my $log = test_log();

my %active_parameter = (
    cluster_constant_path  => catfile(qw{ constant path }),
    conda_path             => catdir( $Bin, qw{ data modules miniconda } ),
    case_id                => 1,
    human_genome_reference => catfile(qw{ a test GRCh37_human_genom_reference.fasta }),
    module_source_environment_command => {
        gatk                   => [ qw{ source activate test_env }, ],
        varianteffectpredictor => [ qw{ source activate test_env }, ],
    },
    outdata_dir                        => catfile(qw{ a outdata dir }),
    program_source_environment_command => {
        picardtools => [ qw{ source activate test_env_1 }, ],
    },
    sample_ids => [qw{ sample_1 }],
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
my @custom_default_parameters =
  qw{ analysis_type bwa_build_reference expansionhunter_repeat_specs_dir exome_target_bed infile_dirs rtg_vcfeval_reference_genome sample_info_file };

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

## Then the defaults should be set for each parameter given
is( $active_parameter{analysis_type}{sample_1}, q{wgs}, q{Set analysis_type default} );

is(
    $active_parameter{bwa_build_reference},
    $active_parameter{human_genome_reference},
    q{Set human_genome_reference default for bwa}
);

is(
    $active_parameter{rtg_vcfeval_reference_genome},
    $active_parameter{human_genome_reference},
    q{Set human_genome_reference default for rtg vcfeval reference genome}
);

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

ok( $active_parameter{expansionhunter_repeat_specs_dir},
    q{Set expansionhunter_repeat_specs_dir} );

## Test setting custom paths
my %test_hash = (
    gatk_path        => catdir( $Bin, qw{ data modules GenomeAnalysisTK-3.7 } ),
    picardtools_path => catdir(
        $active_parameter{conda_path}, qw{ envs test_env_1 share picard-2.14.1-0 }
    ),
    snpeff_path => catdir( $active_parameter{conda_path}, qw{ share snpeff } ),
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
