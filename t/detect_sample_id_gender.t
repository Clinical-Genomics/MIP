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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.03;

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
        q{MIP::Pedigree}       => [qw{ detect_sample_id_gender }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Pedigree qw{ detect_sample_id_gender };

diag(   q{Test detect_sample_id_gender from Pedigree.pm v}
      . $MIP::Pedigree::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given sample ids and genders
my %active_parameter = ( sample_ids => [qw{ sample_1 sample_2 sample_3 }], );

my %sample_info = (
    sample => {
        sample_1 => { sex => q{male}, },
        sample_2 => { sex => q{female}, },
        sample_3 => { sex => q{other}, },
    },
);
(

    $active_parameter{found_male},
    $active_parameter{found_female},
    $active_parameter{found_other},
  )
  = detect_sample_id_gender(
    {
        active_parameter_href => \%active_parameter,
        sample_info_href      => \%sample_info,
    }
  );

my %expected_result = (
    found_male   => 2,
    found_female => 1,
    found_other  => 1,
);

my %expected_gender_info = (
    gender => {
        males   => [qw{ sample_1 }],
        females => [qw{ sample_2 }],
        others  => [qw{ sample_3 }],
    }
);

GENDER:
foreach my $found_gender ( keys %expected_result ) {

## Then all genders count should be one
    is( $active_parameter{$found_gender}, $expected_result{$found_gender},
        $found_gender );
}

## Then sample_ids should be added to each gender category
is_deeply(
    $active_parameter{gender},
    $expected_gender_info{gender},
    q{Added gender info to active parameter}
);

## Given no males or females

%sample_info = (
    sample => {
        sample_1 => { sex => q{xyz}, },
        sample_2 => { sex => q{xyz}, },
        sample_3 => { sex => q{other}, },
    },
);

%expected_result = (
    found_male   => 3,
    found_female => 0,
    found_other  => 3,
);
(

    $active_parameter{found_male},
    $active_parameter{found_female},
    $active_parameter{found_other},
  )
  = detect_sample_id_gender(
    {
        active_parameter_href => \%active_parameter,
        sample_info_href      => \%sample_info,
    }
  );

GENDER:
foreach my $found_gender ( keys %expected_result ) {

## Then one male should be found and a other count of three
    is( $active_parameter{$found_gender}, $expected_result{$found_gender},
        $found_gender );
}

done_testing();
