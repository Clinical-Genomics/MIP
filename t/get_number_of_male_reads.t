#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Path qw{ remove_tree };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp;
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
use MIP::Test::Fixtures qw{ test_log };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Recipes::Analysis::Estimate_gender} => [qw{ get_number_of_male_reads }],
        q{MIP::Test::Fixtures}                     => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Analysis::Estimate_gender qw{ get_number_of_male_reads };

diag(   q{Test get_number_of_male_reads from Estimate_gender.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

test_log( { no_screen => 1, } );

## Given an arbitrary command
my @commands = qw{ ls };

## Given a outscript dir
my $outscript_dir = catfile( File::Temp->newdir() );

## Then return stdout
my $is_ok = get_number_of_male_reads(
    {
        commands_ref  => \@commands,
        outscript_dir => $outscript_dir,
    }
);

ok( $is_ok, q{Executed bash script and captured return} );

remove_tree($outscript_dir);

done_testing();
