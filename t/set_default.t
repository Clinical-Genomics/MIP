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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.01;

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
        q{MIP::Parameter}          => [qw{ set_default }],
        q{MIP::Test::Fixtures}     => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File::Format::Yaml qw{ load_yaml };
use MIP::Parameter qw{ set_default };

diag(   q{Test set_default from Parameter.pm v}
      . $MIP::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Creates log object
my $log = test_log( { no_screen => 1 } );

## Given
my %active_parameter = (
    cluster_constant_path  => catfile(qw{ constant path }),
    conda_path             => catdir( $Bin, qw{ data modules miniconda } ),
    case_id                => 1,
    expansionhunter        => 1,
    frequency_annotation   => 0,
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
    outdata_dir               => catfile(qw{ a outdata dir }),
    rankvariant               => 0,
    sambamba_depth            => 0,
    sample_ids                => [qw{ sample_1 }],
    select_programs           => undef,
    sv_annotate               => 0,
    sv_rankvariant            => 0,
    sv_varianteffectpredictor => 0,
    varianteffectpredictor    => 0,
);

## Mip analyse parameters
## The order of files in @definition_files should follow commands inheritance
my @definition_files =
  ( catfile( dirname($Bin), qw{ definitions analyse_parameters.yaml } ), );

my %parameter;

$parameter{cache}{consensus_analysis_type} = q{wgs};

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

my $is_ok = set_default(
    {
        active_parameter_href => \%active_parameter,
        custom_default_parameters_ref =>
          \@{ $parameter{custom_default_parameters}{default} },
        parameter_href => \%parameter,
    }
);

## Then
ok( $is_ok, q{Set defaults} );

done_testing();
