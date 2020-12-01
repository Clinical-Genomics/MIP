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
use Test::Trap qw{ :stderr:output(systemsafe) };

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_constants test_log test_mip_hashes test_standard_cli };

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
        q{MIP::Analysis}       => [qw{ check_ids_in_dna_vcf }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Analysis qw{ check_ids_in_dna_vcf };

diag(   q{Test check_ids_in_dna_vcf from Analysis.pm v}
      . $MIP::Analysis::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

test_log( {} );

## Given matching dna and rna sample ids
my %active_parameter = test_mip_hashes( { mip_hash_name => q{active_parameter} } );
my %sample_info      = test_mip_hashes( { mip_hash_name => q{qc_sample_info} } );
my $vcf_file_path    = catfile( dirname($Bin),
    qw{ t data test_data 643594-miptest_sorted_md_brecal_comb_BOTH.bcf } );
my %test_process_return = (
    buffers_ref   => [],
    error_message => undef,
    stderrs_ref   => [],
    stdouts_ref   => [q{ADM1059A1 ADM1059A2 ADM1059A3}],
    success       => 1,
);
test_constants(
    {
        test_process_return_href => \%test_process_return,
    }
);

my $is_ok = check_ids_in_dna_vcf(
    {
        active_parameter_href => \%active_parameter,
        dna_vcf_file          => $vcf_file_path,
        sample_info_href      => \%sample_info,
    }
);
## Then return ok
ok( $is_ok, q{Found matching RNA and DNA samples} );

## Given dna samples that doesn't match the rna samples
@{ $active_parameter{sample_ids} } = qw{ ADM1059A4 ADM1059A5 ADM1059A6 };
SAMPLE_ID:
foreach my $sample_id ( @{ $active_parameter{sample_ids} } ) {

    $sample_info{sample}{$sample_id}{subject}{dna_sample_id} = $sample_id;
}
trap {
    check_ids_in_dna_vcf(
        {
            active_parameter_href => \%active_parameter,
            dna_vcf_file          => $vcf_file_path,
            sample_info_href      => \%sample_info,
        }
    )
};

## Then exit and print fatal message
is( $trap->exit, 1, q{Exit when no samples match} );
like(
    $trap->stderr,
    qr/No\smatching\ssample\sids/xms,
    q{Print error message when no samples match}
);

## Given that some dna samples match
@{ $active_parameter{sample_ids} } = qw{ ADM1059A1 ADM1059A5 ADM1059A6 };
$sample_info{sample}{ADM1059A1}{subject}{dna_sample_id} = q{ADM1059A1};
trap {
    check_ids_in_dna_vcf(
        {
            active_parameter_href => \%active_parameter,
            dna_vcf_file          => $vcf_file_path,
            sample_info_href      => \%sample_info,
        }
    )
};
## Then exit and print fatal message
is( $trap->exit, 1, q{Exit on partial sample match} );
like( $trap->stderr, qr/Only\spartial\smatch/xms,
    q{Print error message on partial sample match} );

## Given partial match and force flag
$active_parameter{force_dna_ase} = 1;
trap {
    check_ids_in_dna_vcf(
        {
            active_parameter_href => \%active_parameter,
            dna_vcf_file          => $vcf_file_path,
            sample_info_href      => \%sample_info,
        }
    )
};

## Then warn and set no_ase_samples array
is( $trap->leaveby, q{return}, q{Return on partial sample match and force flag} );
like( $trap->stderr, qr/Turning\soff\sASE/xms,
    q{Print warning message on partial sample match and force flag} );
is_deeply( \@{ $active_parameter{no_ase_samples} },
    [qw{ADM1059A5 ADM1059A6}], q{Set no ASE samples} );

done_testing();
