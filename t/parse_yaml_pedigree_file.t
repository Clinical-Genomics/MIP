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
        q{MIP::Pedigree}       => [qw{ parse_yaml_pedigree_file }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Pedigree qw{ parse_yaml_pedigree_file };

diag(   q{Test parse_yaml_pedigree_file from Pedigree.pm v}
      . $MIP::Pedigree::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

## Given correct input
my %active_parameter = test_mip_hashes( { mip_hash_name => q{active_parameter}, } );
my %parameter        = test_mip_hashes( { mip_hash_name => q{define_parameter}, } );
my %pedigree         = test_mip_hashes( { mip_hash_name => q{pedigree}, } );
my %sample_info;

# Set pedigree file path
$active_parameter{pedigree_file} = catfile( $Bin, qw{ data test_data pedigree.yaml} );

my $is_ok = parse_yaml_pedigree_file(
    {
        active_parameter_href => \%active_parameter,
        file_path             => $active_parameter{pedigree_file},
        parameter_href        => \%parameter,
        pedigree_href         => \%pedigree,
        sample_info_href      => \%sample_info,
    }
);

## Then return true
ok( $is_ok, q{Passed parsing} );

## Given a non matching case_id between cmd and pedigree file
$pedigree{case} = q{not_a_matching_case_id};

trap {
    parse_yaml_pedigree_file(
        {
            active_parameter_href => \%active_parameter,
            file_path             => $active_parameter{pedigree_file},
            parameter_href        => \%parameter,
            pedigree_href         => \%pedigree,
            sample_info_href      => \%sample_info,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if case_id does not match for pedigree and cmd} );
like(
    $trap->stderr,
    qr/for \s+ pedigree \s+ case_id/xms,
    q{Throw fatal log message if case_id does not match for pedigree and cmd}
);

done_testing();
