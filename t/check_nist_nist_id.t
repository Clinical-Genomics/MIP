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
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE};
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

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
        q{MIP::Reference}      => [qw{ check_nist_nist_id }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Reference qw{ check_nist_nist_id };

diag(   q{Test check_nist_nist_id from Reference.pm v}
      . $MIP::Reference::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( {} );

## Given nist info
my %active_parameter = (
    nist_call_set_vcf =>
      { q{3.3.2} => { NA12878 => q{grch37_nist_hg001_-na12878_v3.3.2-.vcf}, }, },
    nist_call_set_bed =>
      { q{3.3.2} => { NA12878 => q{grch37_nist_hg001_-na12878_v3.3.2-.bed}, }, },
    nist_id       => { sample_1 => q{NA12878}, },
    nist_versions => [qw{ 3.3.2 }],
    reference_dir => catdir( $Bin, qw{ data references } ),
    sample_ids    => [qw{ sample_1 }],
);
my @nist_parameters = (qw{ nist_call_set_vcf nist_call_set_bed });

my $is_ok = check_nist_nist_id(
    {
        active_parameter_href => \%active_parameter,
        nist_id_href          => $active_parameter{nist_id},
        nist_parameters_ref   => \@nist_parameters,
    }
);

## Then return true
ok( $is_ok, q{Checked nist nist_id} );

## Given a missing nist_id
$active_parameter{nist_id}{sample_1} = q{not_a_nist_id};

trap {
    check_nist_nist_id(
        {
            active_parameter_href => \%active_parameter,
            nist_id_href          => $active_parameter{nist_id},
            nist_parameters_ref   => \@nist_parameters,
        }
    )
};

## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if nist_id is not mapped} );
like( $trap->stderr, qr/Supplied\s+nist\s+id:\s+/xms, q{Throw fatal log message} );

done_testing();
