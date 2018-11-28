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
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

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
            sample_origin     => q{Valhalla},
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
    sample_origin     => q{Valhalla},
    sex               => q{male},
);

while ( my ( $attribute, $attribute_value ) = each %attribute ) {

## Then
    my $got_attribute = get_pedigree_sample_id_attributes(
        {
            attribute        => $attribute,
            sample_id        => $sample_id,
            sample_info_href => \%sample_info,
        }
    );

    is( $got_attribute, $attribute_value,
            q{Got }
          . $sample_id
          . q{ attribute: }
          . $attribute . q{ => }
          . $attribute_value );
}

done_testing();
