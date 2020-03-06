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
        q{MIP::Io::Recipes}    => [qw{ set_file_prefix_tag }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Io::Recipes qw{ set_file_prefix_tag };

diag(   q{Test set_file_prefix_tag from Recipes.pm v}
      . $MIP::Io::Recipes::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $sample_id     = q{sample-1};
my $current_chain = q{MAIN};

my @order_parameters = qw{ bwa_mem pmerge pmark };
my %active_parameter = (
    bwa_mem => 1,
    manta   => 1,
    pmark   => 1,
    pmerge  => 0,
);
my %expected_file_tag;
my %file_info;
my %file_tag = (
    bwa_mem => q{mem},
    manta   => q{manta},
    pmark   => q{md},
    pmerge  => q{merge},
);
my %temp_file_ending;

## Given MAIN chain file tags, when merge is turned off
RECIPE:
foreach my $recipe_name (@order_parameters) {

    $temp_file_ending{$current_chain}{$sample_id} = set_file_prefix_tag(
        {
            current_chain         => $current_chain,
            file_tag              => $file_tag{$recipe_name},
            file_info_href        => \%file_info,
            id                    => $sample_id,
            is_active_recipe      => $active_parameter{$recipe_name},
            recipe_name           => $recipe_name,
            temp_file_ending_href => \%temp_file_ending,
        }
    );

}

## Define what to expect
# First file tag
$expected_file_tag{$sample_id}{bwa_mem}{file_tag} = $file_tag{bwa_mem};

# Propagated
$expected_file_tag{$sample_id}{pmerge}{file_tag} = $file_tag{bwa_mem};

# Sequential
$expected_file_tag{$sample_id}{pmark}{file_tag} =
  $file_tag{bwa_mem} . $file_tag{pmark};

## Then 3 file tags should be added where one is sequential and one is just propagated
is_deeply( \%file_info, \%expected_file_tag, q{Added file prefix tags for MAIN } );

## Given other chain than MAIN
push @order_parameters, q{manta};
my $other_chain = q{SV};

# Clear previous file tag builds
%file_info        = ();
%temp_file_ending = ();

## Given SV chain file tags, when merge is turned off
RECIPE:
foreach my $recipe_name (@order_parameters) {

    $temp_file_ending{$current_chain}{$sample_id} = set_file_prefix_tag(
        {
            current_chain         => $current_chain,
            file_tag              => $file_tag{$recipe_name},
            file_info_href        => \%file_info,
            id                    => $sample_id,
            is_active_recipe      => $active_parameter{$recipe_name},
            recipe_name           => $recipe_name,
            temp_file_ending_href => \%temp_file_ending,
        }
    );
}
$expected_file_tag{$sample_id}{manta}{file_tag} =
  $expected_file_tag{$sample_id}{pmark}{file_tag} . $file_tag{manta};

is_deeply( \%file_info, \%expected_file_tag,
    q{Added file prefix tags for chain that inherits from MAIN } );

done_testing();
