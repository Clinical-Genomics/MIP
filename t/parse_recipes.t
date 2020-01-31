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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

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
        q{MIP::Recipes::Parse} => [qw{ parse_recipes }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Parse qw{ parse_recipes };

diag(   q{Test parse_recipes from Parse.pm v}
      . $MIP::Recipes::Parse::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { no_screen => 0 } );

## Given recipes hashes
my %active_parameter = ( recipe_time => { bwa_mem => 1, }, );

my %parameter = (
    bwa_mem => {
        analysis_mode       => q{sample},
        associated_recipe   => [qw{ mip }],
        data_type           => q{SCALAR},
        default             => 1,
        file_tag            => q{_sorted},
        outfile_suffix      => q{.bam},
        program_executables => [qw{ bwa samtools sambamba }],
        recipe_type         => q{aligners},
        type                => q{recipe},
    },
    mip => {
        associated_recipe => [ qw{ mip }, ],
        data_type         => q{SCALAR},
        default           => 1,
        type              => q{mip}
    },
    recipe_time => {
        associated_recipe => [ qw{ mip }, ],
        data_type         => q{HASH},
        default           => { bwa_mem => 1, },
        type              => q{mip},
    },
);
my %parameter_to_check = (
    keys     => [qw{ recipe_time  }],
    elements => [qw{ associated_recipe decompose_normalize_references }],
);

my $is_ok = parse_recipes(
    {
        active_parameter_href   => \%active_parameter,
        parameter_href          => \%parameter,
        parameter_to_check_href => \%parameter_to_check,
    }
);

## Then return true if parsing passed
ok( $is_ok, q{Parsed recipe names} );

done_testing();
