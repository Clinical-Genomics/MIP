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
use Readonly;
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $EMPTY_STR $SPACE };
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
        q{MIP::Analysis}       => [qw{ update_recipe_mode_for_analysis_type }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Analysis qw{ update_recipe_mode_for_analysis_type };

diag(   q{Test update_recipe_mode_for_analysis_type from Analysis.pm v}
      . $MIP::Analysis::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Create log object
my $log = test_log( {} );

my @recipes = qw{ cnvnator_ar delly_call delly_reformat samtools_subsample_mt tiddit };
my %active_parameter = (
    cnvnator_ar           => 1,
    delly_call            => 1,
    delly_reformat        => 1,
    manta                 => 1,
    samtools_subsample_mt => 1,
    tiddit                => 1,
);

trap {
    update_recipe_mode_for_analysis_type(
        {
            active_parameter_href   => \%active_parameter,
            consensus_analysis_type => q{wgs},
            recipes_ref             => \@recipes,
        }
    )
};

is( $trap->stderr, $EMPTY_STR, q{No updates to recipes mode} );

trap {
    update_recipe_mode_for_analysis_type(
        {
            active_parameter_href   => \%active_parameter,
            consensus_analysis_type => q{wes},
            recipes_ref             => \@recipes,
        }
    )
};
## Unpack
my $cnvnator_mode              = $active_parameter{cnvnator_ar};
my $delly_call_mode            = $active_parameter{delly_call};
my $delly_reformat_mode        = $active_parameter{delly_reformat};
my $manta_mode                 = $active_parameter{manta};
my $samtools_subsample_mt_mode = $active_parameter{samtools_subsample_mt};
my $tiddit_mode                = $active_parameter{tiddit};

## Test recipe mode updates and warnings
is( $cnvnator_mode,              0, q{Updated recipes mode for cnvnator_ar} );
is( $delly_call_mode,            0, q{Updated recipes mode for delly_call} );
is( $delly_reformat_mode,        0, q{Updated recipes mode for delly_reformat} );
is( $manta_mode,                 1, q{Updated recipes mode for manta} );
is( $samtools_subsample_mt_mode, 0, q{Updated recipes mode for samtools_subsample_mt} );
is( $tiddit_mode,                0, q{Updated recipes mode for tiddit} );
like( $trap->stderr, qr/WARN/xms, q{Generated warning message} );

done_testing();
