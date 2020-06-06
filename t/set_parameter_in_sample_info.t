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
use MIP::Constants qw{ $COMMA $NEWLINE $SPACE };
use MIP::Test::Fixtures qw{ test_mip_hashes test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.02;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $GENOME_BUILD_VERSION => 37;

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Sample_info}    => [qw{ set_parameter_in_sample_info }],
        q{MIP::Test::Fixtures} => [qw{ test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Sample_info qw{ set_parameter_in_sample_info };

diag(   q{Test set_parameter_in_sample_info from Sample_info.pm v}
      . $MIP::Sample_info::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given parameter paths
my %active_parameter = test_mip_hashes(
    {
        mip_hash_name => q{active_parameter},
    }
);
%{ $active_parameter{expected_coverage} } = (
    ADM1059A1 => 1,
    ADM1059A2 => 1,
    ADM1059A3 => 1,
);
$active_parameter{human_genome_reference} = catfile(qw{ a test path genome build});
$active_parameter{log_file}               = catfile(qw{ a test dir and log_path});
$active_parameter{pedigree_file}          = catfile(qw{ a test pedigree path });

my %sample_info;

my %file_info = (
    human_genome_reference_source  => q{grch},
    human_genome_reference_version => $GENOME_BUILD_VERSION,
);

set_parameter_in_sample_info(
    {
        active_parameter_href => \%active_parameter,
        file_info_href        => \%file_info,
        sample_info_href      => \%sample_info,
    }
);

my %expected_sample_info = (
    human_genome_build => {
        path    => $active_parameter{human_genome_reference},
        source  => $file_info{human_genome_reference_source},
        version => $file_info{human_genome_reference_version},
    },
);
$expected_sample_info{analysis_type} = $active_parameter{analysis_type};
%{ $expected_sample_info{expected_coverage} } = %{ $active_parameter{expected_coverage} };
$expected_sample_info{pedigree_file}{path} = $active_parameter{pedigree_file};
$expected_sample_info{log_file_dir} = dirname( dirname( $active_parameter{log_file} ) );
$expected_sample_info{last_log_file_path} = $active_parameter{log_file};
$expected_sample_info{has_trio}           = 0;

## Then these entries should be set in sample info
is_deeply( \%sample_info, \%expected_sample_info, q{Set parameter in sample_info} );

done_testing();
