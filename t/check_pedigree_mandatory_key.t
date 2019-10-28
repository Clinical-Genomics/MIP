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
use Modern::Perl qw{ 2018 };
use Readonly;
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

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
        q{MIP::Check::Pedigree} => [qw{ check_pedigree_mandatory_key }],
        q{MIP::Test::Fixtures}  => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Check::Pedigree qw{ check_pedigree_mandatory_key };

diag(   q{Test check_pedigree_mandatory_key from Pedigree.pm v}
      . $MIP::Check::Pedigree::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Creates log object
my $log = tets_log( {} );

my $active_parameter{dna_vcf_file} = catfile(q{dna.vcf});

my %pedigree = (
    case    => q{case_1},
    samples => [
        {
            analysis_type => q{wes},
            dna_sample_id => q{sample_1},
            father        => 0,
            mother        => 0,
            phenotype     => q{affected},
            sample_id     => q{sample_1},
            sex           => q{female},
        },
        {
            analysis_type => q{wgs},
            dna_sample_id => q{sample_2},
            father        => 0,
            mother        => 0,
            phenotype     => q{unaffected},
            sample_id     => q{sample_2},
            sex           => q{male},
        },
        {
            analysis_type => q{wts},
            dna_sample_id => q{sample_3},
            father        => 0,
            mother        => 0,
            phenotype     => q{unknown},
            sample_id     => q{sample_3},
            sex           => q{other},
        },
        {
            analysis_type => q{wgs},
            dna_sample_id => q{sample_4},
            father        => q{sample_1},
            mother        => q{sample_2},
            phenotype     => q{unknown},
            sample_id     => q{sample_4},
            sex           => q{unknown},
        },
    ],
);
##Given all mandatory keys
my $is_ok = check_pedigree_mandatory_key(
    {
        active_parameter_href => \%active_parameter,
        file_path             => catfile(qw{ path to pedigree.yaml }),
        pedigree_href         => \%pedigree,
    }
);

## Then return true
ok( $is_ok, q{Found all mandatory keys in pedigree yaml file} );

## Given a missing mandatory case key
delete $pedigree{case};

trap {
    check_pedigree_mandatory_key(
        {
            active_parameter_href => \%active_parameter,
            file_path             => catfile(qw{ path to pedigree.yaml }),
            pedigree_href         => \%pedigree,
        }
      )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if mandatory case key cannot be found} );
like( $trap->stderr, qr/FATAL/xms,
    q{Throw fatal log message if mandatory case key cannot be found} );

# Reinitialize key
$pedigree{case} = q{case_1};

## Given a missing mandatory sample key
delete $pedigree{samples}[0]{phenotype};

trap {
    check_pedigree_mandatory_key(
        {
            active_parameter_href => \%active_parameter,
            file_path             => catfile(qw{ path to pedigree.yaml }),
            pedigree_href         => \%pedigree,
        }
      )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if mandatory sample key cannot be found} );
like( $trap->stderr, qr/FATAL/xms,
    q{Throw fatal log message if mandatory sample key cannot be found} );

done_testing();
