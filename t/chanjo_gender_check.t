#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir };
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
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_mip_hashes test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

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
        q{MIP::Check::Qccollect} => [qw{ chanjo_gender_check }],
        q{MIP::Test::Fixtures}   => [qw{ test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Check::Qccollect qw{ chanjo_gender_check };

diag(   q{Test chanjo_gender_check from Qccollect.pm v}
      . $MIP::Check::Qccollect::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given pedigree and sample_info sex, when represented by letters
my $chanjo_sexcheck_gender = q{female};
my %expected_qc_data;
my $infile = q{an_infile};
my %qc_data;
my %sample_info = test_mip_hashes(
    {
        mip_hash_name => q{qc_sample_info},
        recipe_name   => q{a_recipe},
    }
);
my %pedigree = (
    ADM1059A1 => q{female},
    ADM1059A2 => q{male},
    ADM1059A3 => q{unknown},
);

while ( my ( $sample_id, $gender ) = each %pedigree ) {

    chanjo_gender_check(
        {
            chanjo_sexcheck_gender => $gender,
            infile                 => $infile,
            qc_data_href           => \%qc_data,
            sample_id              => $sample_id,
            sample_info_href       => \%sample_info,
        }
    );
    $expected_qc_data{sample}{$sample_id}{$infile}{gender_check} = q{PASS};

## Then PASS should be added to qc_data for sample
    is_deeply( \%qc_data, \%expected_qc_data, q{Passed } . $gender . q{ letters} );
}

## Given pedigree and sample_info sex, when represented by numbers
# Modify to plink format for gender
$sample_info{sample}{ADM1059A1}{sex} = 2;
$sample_info{sample}{ADM1059A2}{sex} = 1;
$sample_info{sample}{ADM1059A3}{sex} = 0;

while ( my ( $sample, $gender ) = each %pedigree ) {

    chanjo_gender_check(
        {
            chanjo_sexcheck_gender => $gender,
            infile                 => $infile,
            qc_data_href           => \%qc_data,
            sample_id              => $sample,
            sample_info_href       => \%sample_info,
        }
    );
    $expected_qc_data{sample}{$sample}{$infile}{gender_check} = q{PASS};

## Then PASS should be added to qc_data for sample
    is_deeply( \%qc_data, \%expected_qc_data, q{Passed } . $gender . q{ numbers} );
}

## Given an incorrect gender
my $failed_sample_id = q{ADM1059A1};

chanjo_gender_check(
    {
        chanjo_sexcheck_gender => q{male},
        infile                 => $infile,
        qc_data_href           => \%qc_data,
        sample_id              => $failed_sample_id,
        sample_info_href       => \%sample_info,
    }
);
$expected_qc_data{sample}{$failed_sample_id}{$infile}{gender_check} = q{FAIL};

is_deeply( \%qc_data, \%expected_qc_data, q{Failed gender } );

done_testing();
