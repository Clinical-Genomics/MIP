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
use MIP::Test::Fixtures qw{ test_constants test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.02;

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
        q{MIP::Vcf}            => [qw{ get_sample_ids_from_vcf }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Vcf qw{ get_sample_ids_from_vcf };

diag(   q{Test get_sample_ids_from_vcf from Vcf.pm v}
      . $MIP::Vcf::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

test_log( {} );
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

## Given vcf input
my $vcf_file_path = catfile( dirname($Bin),
    qw{ t data test_data 643594-miptest_sorted_md_brecal_comb_BOTH.bcf } );
my @sample_ids = get_sample_ids_from_vcf(
    {
        vcf_file_path => $vcf_file_path,
    }
);

## Then return sample ids
my @expected_ids = qw{ ADM1059A1 ADM1059A2 ADM1059A3 };
is_deeply( \@sample_ids, \@expected_ids, q{Get VCF sample ids} );

## Given something that will error
$test_process_return{stderrs_ref} = [q{error}];
test_constants(
    {
        test_process_return_href => \%test_process_return,
    }
);
trap {
    get_sample_ids_from_vcf(
        {
            vcf_file_path => catfile(qw{ not a file.vcf }),
        }
    )
};

## Then exit and print fatal message
is( $trap->exit, 1, q{Exit on error} );
like( $trap->stderr, qr/ERROR:/xms, q{Print error message} );

done_testing();
