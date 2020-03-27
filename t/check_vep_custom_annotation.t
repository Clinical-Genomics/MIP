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
        q{MIP::Vep}            => [qw{ check_vep_custom_annotation }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Vep qw{ check_vep_custom_annotation };

diag(   q{Test check_vep_custom_annotation from Vep.pm v}
      . $MIP::Vep::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

## Given a undefined vep_custom_annotation
my %vep_custom_ann_undef;
my $is_not_ok = check_vep_custom_annotation(
    {
        vep_custom_ann_href => \%vep_custom_ann_undef,
    }
);

## Then there is nothing to check - return false
is( $is_not_ok, 0, q{Nothing to check} );

## Given a valid annotation
my %vep_custom_ann = (
    genomic_superdups_frac_match => {
        annotation_type          => q{overlap},
        force_report_coordinates => 0,
        key                      => q{genomic_superdups_frac_match},
        file_type                => q{bed},
        path                     => catfile(
            $Bin, qw{ data references grch37_genomics_super_dups_-20181009.bed.gz }
        ),
    },
);

my $is_ok = check_vep_custom_annotation(
    {
        vep_custom_ann_href => \%vep_custom_ann,
    }
);
## Then return true
ok( $is_ok, q{Passed checks} );

## Given a not valid hash ref
$vep_custom_ann{not_valid_annotation} = [q{not_a_valid_ref}];

trap {
    check_vep_custom_annotation(
        {
            vep_custom_ann_href => \%vep_custom_ann,
        }
    )
};

## Then exit and throw FATAL log message
is( $trap->leaveby, q{die}, q{Exit if not a hash ref } );
like( $trap->die, qr/Is\s+not\s+a/xms, q{Not a hash ref} );

delete $vep_custom_ann{not_valid_annotation};

## Given not valid annotation options when no path supplied
my %vep_custom_ann_not_valid_option = (
    genomic_superdups_frac_match => {
        annotation_type          => q{not_an_annotation_type},
        force_report_coordinates => q{not_a_boolean},
        file_type                => q{not_allowed_file_type},
        key                      => q{genomic_superdups_frac_match},
    },
);
trap {
    check_vep_custom_annotation(
        {
            vep_custom_ann_href => \%vep_custom_ann_not_valid_option,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if not valid options value for annotation } );
like(
    $trap->stderr,
    qr/has\s+a\s+not\s+allowed\s+option/xms,
    q{Not valid options for annotation}
);

## Given valid options except path option
my %vep_custom_ann_bad_path = (
    genomic_superdups_frac_match => {
        annotation_type          => q{overlap},
        force_report_coordinates => 0,
        key                      => q{genomic_superdups_frac_match},
        file_type                => q{bed},
        path                     => catfile( $Bin, qw{ not a path } ),
    },
);
trap {
    check_vep_custom_annotation(
        {
            vep_custom_ann_href => \%vep_custom_ann_bad_path,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if path does not exist for annotation } );
like( $trap->stderr, qr/FATAL/xms, q{Throw fatal log message if no path} );

## Given a missing key option
my %vep_custom_ann_bad_key = (
    genomic_superdups_frac_match => {
        annotation_type          => q{overlap},
        force_report_coordinates => 0,
        file_type                => q{bed},
        path                     => catfile( $Bin, qw{ not a path } ),
    },
);
trap {
    check_vep_custom_annotation(
        {
            vep_custom_ann_href => \%vep_custom_ann_bad_key,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if key does not exist for annotation } );
like(
    $trap->stderr,
    qr/lacks\s+required\s+option/xms,
    q{Throw fatal log message if no key for annotation}
);

done_testing();
