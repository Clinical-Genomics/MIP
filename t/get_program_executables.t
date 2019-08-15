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
use MIP::Test::Fixtures qw{ test_mip_hashes test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Get::Parameter} => [qw{ get_program_executables }],
        q{MIP::Test::Fixtures} => [qw{ test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Get::Parameter qw{ get_program_executables };

diag(   q{Test get_program_executables from Parameter.pm v}
      . $MIP::Get::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given
my %parameter = test_mip_hashes( { mip_hash_name => q{define_parameter}, } );

trap {
    get_program_executables( { parameter_href => \%parameter, } )
};

## Then all program executables is returned
## Then croak and exist
is( $trap->leaveby, q{die}, q{Exit if no recipe is found} );
like(
    $trap->die,
    qr{ No\skeys\s'cache'\sand\s'recipe'}xms,
    q{Throw error if no recipe can be found}
);

## Given a cache recipe and program
push @{ $parameter{cache}{recipe} }, q{bwa_mem};

my @program_executables = get_program_executables( { parameter_href => \%parameter, } );

my @expected_program_executables = qw{ bwa samtools sambamba };

## Then all program executables is returned
is_deeply( \@expected_program_executables,
    \@program_executables, q{Got all program executables} );

done_testing();
