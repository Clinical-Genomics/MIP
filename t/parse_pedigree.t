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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_mip_hashes };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Pedigree}       => [qw{ parse_pedigree }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_mip_hashes }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Pedigree qw{ parse_pedigree };

diag(   q{Test parse_pedigree from Pedigree.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { no_screen => 1, } );

## Given correct input
my %active_parameter = test_mip_hashes( { mip_hash_name => q{active_parameter}, } );
my %parameter        = test_mip_hashes( { mip_hash_name => q{define_parameter}, } );
my %sample_info;

# Set pedigree file path
$active_parameter{pedigree_file} = catfile( $Bin, qw{ data test_data pedigree_wes.yaml} );

my $is_ok = parse_pedigree(
    {
        active_parameter_href => \%active_parameter,
        pedigree_file_path    => $active_parameter{pedigree_file},
        parameter_href        => \%parameter,
        sample_info_href      => \%sample_info,
    }
);

## Then return true
ok( $is_ok, q{Passed parsing of pedigree data} );

done_testing();
