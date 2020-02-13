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
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

my $VERBOSE = 1;
our $VERSION = '1.0.0';

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $COLON => q{:};
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module =
      ( q{MIP::Test::Fixtures} => [qw{ test_log test_mip_hashes test_standard_cli }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Test::Fixtures qw{ test_mip_hashes };

diag(   q{Test test_mip_hashes from Fixtures.pm v}
      . $MIP::Test::Fixtures::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

## Given active parameters
my $recipe_name = q{bwa_mem};

my %active_parameter = test_mip_hashes(
    {
        mip_hash_name => q{active_parameter},
        recipe_name   => $recipe_name,
    }
);

## Then dynamic parameters should be set
is( $active_parameter{$recipe_name}, 2, q{Set recipe mode} );
ok( $active_parameter{temp_directory}, q{Set temp_directory} );

## Given file info parameters
my %file_info = test_mip_hashes(
    {
        mip_hash_name => q{file_info},
        recipe_name   => $recipe_name,
    }
);
## Then dynamic parameters should be set
is(
    $file_info{human_genome_reference},
    q{grch37_homo_sapiens_-d5-.fasta},
    q{Set human genome reference}
);

## Given parameters
my %parameter = test_mip_hashes( { mip_hash_name => q{recipe_parameter}, } );

## Then dynamic parameters should be set
is( $parameter{$recipe_name}{chain}, q{TEST}, q{Set recipe chain} );

done_testing();
