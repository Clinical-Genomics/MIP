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
        q{MIP::Check::Parameter} => [qw{ check_recipe_fastq_compatibility }],
        q{MIP::Test::Fixtures}   => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Check::Parameter qw{ check_recipe_fastq_compatibility };
use MIP::Test::Fixtures qw{ test_log test_mip_hashes };

diag(   q{Test check_recipe_fastq_compatibility from Parameter.pm v}
      . $MIP::Check::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my %dependency_tree = test_mip_hashes(
    {
        mip_hash_name => q{dependency_tree_rna},
    }
);
my %parameter = test_mip_hashes(
    {
        mip_hash_name => q{define_parameter},
    }
);
my %active_parameter = test_mip_hashes(
    {
        mip_hash_name => q{active_parameter},
    }
);

$parameter{dependency_tree_href} = \%dependency_tree;

my $log = test_log( {} );

my $sample_id          = q{test};
my %infile_lane_prefix = ( $sample_id => [qw{ test_lane1 test_lane2 }], );
my %sample_info        = (
    sample => {
        $sample_id => {
            file => {
                test_lane1 => {
                    sequence_run_type => q{paired-end},
                },
                test_lane2 => {
                    sequence_run_type => q{paired-end},
                },
            },
        },
    },
);

## Given that both lanes have been sequenced the same way
my $recipe_fastq_compatibility = check_recipe_fastq_compatibility(
    {
        active_parameter_href   => \%active_parameter,
        infile_lane_prefix_href => \%infile_lane_prefix,
        parameter_href          => \%parameter,
        recipe_name             => q{salmon_quant},
        sample_info_href        => \%sample_info,
    }
);
## Then OK
ok( $recipe_fastq_compatibility, q{Compatible} );

## Given a lane difference
$sample_info{sample}{$sample_id}{file}{test_lane2}{sequence_run_type} = q{single_end};

trap {
    $recipe_fastq_compatibility = check_recipe_fastq_compatibility(
        {
            active_parameter_href   => \%active_parameter,
            infile_lane_prefix_href => \%infile_lane_prefix,
            parameter_href          => \%parameter,
            recipe_name             => q{salmon_quant},
            sample_info_href        => \%sample_info,
        }
    )
};
## Then not compatible
is( $recipe_fastq_compatibility, 0, q{Identify non compatible sequence types} );
like( $trap->stderr, qr/Multiple\ssequence\srun\stypes\sdetected/xms, q{Log warning} );

done_testing();
