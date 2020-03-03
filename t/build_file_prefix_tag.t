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
        q{MIP::File::Format::Mip} => [qw{ build_file_prefix_tag }],
        q{MIP::Test::Fixtures}    => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File::Format::Mip qw{ build_file_prefix_tag };

diag(   q{Test build_file_prefix_tag from Mip.pm v}
      . $MIP::File::Format::Mip::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $sample_id     = q{homer};
my $case_id       = q{simpsons};
my $current_chain = q{MAIN};
my $other_chain   = q{SV};

my @order_recipes = qw{ bwa_mem pmerge pmark manta };

my %active_parameter = (
    case_id               => $case_id,
    bwa_mem               => 1,
    manta                 => 1,
    pmark                 => 1,
    pmerge                => 0,
    random_test_parameter => undef,
    random_test_parameter => q{not_a_recipe},
    sample_ids            => [$sample_id],
);
my %file_info;
my %parameter = (
    cache   => { recipe => [qw{ bwa_mem pmark manta pmerge }], },
    bwa_mem => {
        chain    => $current_chain,
        file_tag => q{mem},
    },
    pmark => {
        chain    => $current_chain,
        file_tag => q{md},
    },
    manta => {
        chain    => $other_chain,
        file_tag => q{manta},
    },
    pmerge => {
        chain    => $current_chain,
        file_tag => q{nofile_tag},
    },
);

## Given file tags, when multiple chains
build_file_prefix_tag(
    {
        active_parameter_href => \%active_parameter,
        file_info_href        => \%file_info,
        order_recipes_ref     => \@order_recipes,
        parameter_href        => \%parameter,
    }
);

my %expected_file_tag = (
    $sample_id => {
        bwa_mem => { file_tag => q{mem}, },
        pmark   => { file_tag => q{memmd}, },
        manta   => { file_tag => q{memmdmanta}, },
    },
    $case_id => {
        bwa_mem => { file_tag => q{mem}, },
        pmark   => { file_tag => q{memmd}, },
        manta   => { file_tag => q{memmdmanta}, },
    },
);

## Then these file tags should be set according to %expected_file_tag
is_deeply( \%file_info, \%expected_file_tag, q{Built file endings} );

done_testing();
