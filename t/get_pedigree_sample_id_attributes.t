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
use MIP::Test::Fixtures qw{ test_standard_cli };

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
        q{MIP::Get::Parameter} => [qw{ get_pedigree_sample_id_attributes }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Get::Parameter qw{ get_pedigree_sample_id_attributes };

diag(   q{Test get_pedigree_sample_id_attributes from Parameter.pm v}
      . $MIP::Get::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a sample with info
my $sample_id   = q{God of thunder};
my %sample_info = (
    sample => {
        $sample_id => {
            analysis_type     => q{wgs},
            capture_kit       => q{agilent_v1},
            expected_coverage => 1,
            father            => q{Odin},
            mother            => q{Frey},
            phenotype         => q{affected},
            sample_id         => q{God of thunder},
            sample_name       => q{Thor},
            is_from_sample    => q{Valhalla},
            sex               => q{male},
        },
    },
);
my %attribute = (
    analysis_type     => q{wgs},
    capture_kit       => q{agilent_v1},
    expected_coverage => 1,
    father            => q{Odin},
    mother            => q{Frey},
    phenotype         => q{affected},
    sample_id         => q{God of thunder},
    sample_name       => q{Thor},
    is_from_sample    => q{Valhalla},
    sex               => q{male},
);

while ( my ( $attribute, $attribute_value ) = each %attribute ) {

    my $got_attribute = get_pedigree_sample_id_attributes(
        {
            attribute        => $attribute,
            sample_id        => $sample_id,
            sample_info_href => \%sample_info,
        }
    );

    ## Then we should get an attribute
    is( $got_attribute, $attribute_value,
            q{Got }
          . $sample_id
          . q{ attribute: }
          . $attribute . q{ => }
          . $attribute_value );
}

## Given a undefined attibute, when not mandatory
delete $sample_info{sample}{$sample_id}{expected_coverage};

my $got_attribute = get_pedigree_sample_id_attributes(
    {
        attribute        => q{expected_coverage},
        sample_id        => $sample_id,
        sample_info_href => \%sample_info,
    }
);

## Then return false
is( $got_attribute, undef, q{Returned undef for not mandatory pedigree attribute} );

## Given a undefined attibute, when mandatory
delete $sample_info{sample}{$sample_id}{sex};

trap {
    $got_attribute = get_pedigree_sample_id_attributes(
        {
            attribute        => q{sex},
            sample_id        => $sample_id,
            sample_info_href => \%sample_info,
        }
    )
};

## Then exit and throw FATAL log message
is( $trap->leaveby, q{die}, q{Exit if the attibute is undef} );
like( $trap->die, qr/Could\s+not\s+find\s+sample_info_name/xms, q{Throw error} );

done_testing();
