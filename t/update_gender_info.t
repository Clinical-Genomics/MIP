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
use MIP::Test::Fixtures qw{ test_log test_mip_hashes };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Recipes::Analysis::Estimate_gender} => [qw{ update_gender_info }],
        q{MIP::Test::Fixtures}                     => [qw{ test_log test_mip_hashes }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Analysis::Estimate_gender qw{ update_gender_info };

diag(   q{Test update_gender_info from Estimate_gender.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

test_log( { no_screen => 1, } );

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
my $consensus_analysis_type = q{wgs};
my %file_info               = test_mip_hashes(
    {
        mip_hash_name => q{file_info},
    }
);
my %sample_info = test_mip_hashes(
    {
        mip_hash_name => q{qc_sample_info},
    }
);

## Given a sample id
my $sample_id = $active_parameter{sample_ids}[2];

## Given a male read count
my $y_read_count = $MALE_THRESHOLD + 1;

## Given an unknown gender for a sample
$active_parameter{gender}{others} = [$sample_id];

## When updating gender info
my $is_ok = update_gender_info(
    {
        active_parameter_href   => \%active_parameter,
        consensus_analysis_type => $consensus_analysis_type,
        file_info_href          => \%file_info,
        sample_id               => $sample_id,
        sample_info_href        => \%sample_info,
        y_read_count            => $y_read_count,
    }
);

## Then return true
ok( $is_ok, q{Updated gender info} );

## Then set include_y to 1
is( $active_parameter{include_y}, 1, q{Include y} );

## Then add estimated gender for sample id
is( $active_parameter{gender_estimation}{$sample_id}, q{male}, q{Added estimated gender male} );

## Then add the sample_id to list of male samples ids
is_deeply( $active_parameter{gender}{males}, [$sample_id], q{Add to males sample_id} );

## Then remove the sample_id from list of other genders
is( @{ $active_parameter{gender}{others} }, 0, q{Remove from others sample_id} );

## Given a y read count when female
$y_read_count = $FEMALE_THRESHOLD;
delete $active_parameter{gender};

## Given an unknown gender for a sample
$active_parameter{gender}{others} = [$sample_id];

## When updating gender info
update_gender_info(
    {
        active_parameter_href   => \%active_parameter,
        consensus_analysis_type => $consensus_analysis_type,
        file_info_href          => \%file_info,
        sample_id               => $sample_id,
        sample_info_href        => \%sample_info,
        y_read_count            => $y_read_count,
    }
);

## Then set include y to zero
is( $active_parameter{include_y}, 0, q{Exclude y} );

## Then add estimated gender for sample id
is( $active_parameter{gender_estimation}{$sample_id}, q{female}, q{Added estimated gender female} );

## Then add the sample_id to list of female samples ids
is_deeply( $active_parameter{gender}{females}, [$sample_id], q{Add to females sample_id} );

## Then remove the sample_id from list of other genders
is( @{ $active_parameter{gender}{others} }, 0, q{Remove from others sample_id} );

done_testing();
