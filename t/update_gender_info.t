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
use MIP::Test::Fixtures qw{ test_log test_mip_hashes test_standard_cli };

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
        q{MIP::Parse::Gender}  => [qw{ update_gender_info }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Parse::Gender qw{ update_gender_info };

diag(   q{Test update_gender_info from Gender.pm v}
      . $MIP::Parse::Gender::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { no_screen => 1, } );

## Constants
Readonly my $MALE_THRESHOLD   => 36;
Readonly my $FEMALE_THRESHOLD => 24;

## Given a y read count when male
my %active_parameter = test_mip_hashes(
    {
        mip_hash_name => q{active_parameter},
        recipe_name   => q{bwa_mem},
    }
);

my %file_info = test_mip_hashes(
    {
        mip_hash_name => q{file_info},
    }
);
my $sample_id    = $active_parameter{sample_ids}[2];
my $y_read_count = $MALE_THRESHOLD + 1;

my $is_ok = update_gender_info(
    {
        active_parameter_href => \%active_parameter,
        file_info_href        => \%file_info,
        sample_id             => $sample_id,
        y_read_count          => $y_read_count,
    }
);

## Then return true
ok( $is_ok, q{Updated gender info} );

## Then increment found male and add estimated gender for sample id
is( $active_parameter{found_male}, 1, q{Incremented found male} );
is( $active_parameter{gender_estimation}{$sample_id},
    q{male}, q{Added estimated gender male} );

## Given a y read count when female
$y_read_count = $FEMALE_THRESHOLD;

update_gender_info(
    {
        active_parameter_href => \%active_parameter,
        file_info_href        => \%file_info,
        sample_id             => $sample_id,
        y_read_count          => $y_read_count,
    }
);

## Then increment found female and add estimated gender for sample id
is( $active_parameter{found_male}, 0, q{Decremented found male} );
is( $active_parameter{gender_estimation}{$sample_id},
    q{female}, q{Added estimated gender female} );

done_testing();
