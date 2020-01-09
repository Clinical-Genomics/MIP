#!/usr/bin/env perl

use 5.026;
use Carp;
use Cwd qw{ abs_path };
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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli };

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
        q{MIP::Active_parameter} => [qw{ update_to_absolute_path }],
        q{MIP::Test::Fixtures}   => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Active_parameter qw{ update_to_absolute_path };

diag(   q{Test update_to_absolute_path from Active_parameter.pm v}
      . $MIP::Active_parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given different type of parameters to update path for
my %parameter = (
    hash   => { update_path => q{absolute_path}, },
    array  => { update_path => q{absolute_path}, },
    scalar => { update_path => q{absolute_path}, },
);

my %active_parameter = (
    hash => {
        file => q{annotation},
    },
    array => [ catfile( $Bin, qw{ data references grch37_homo_sapiens_-d5-.fasta.gz } ) ],
    scalar => catfile( $Bin, qw{ data references grch37_homo_sapiens_-d5-.fasta.gz } ),
);

update_to_absolute_path(
    {
        parameter_href        => \%parameter,
        active_parameter_href => \%active_parameter,
    }
);

my $expected_absolute_path =
  catfile( $Bin, qw{ data references grch37_homo_sapiens_-d5-.fasta.gz } );

my $expeceted_hash_key_absolute_path = abs_path( catfile(q{file}) );

## Then update_to_absolute_path uppdated path for hash key and not value
foreach my $key ( keys %{ $active_parameter{hash} } ) {

    is( $key, $expeceted_hash_key_absolute_path, q{Updated hash key to absolute path} );
}

## Then elements path should now be absolute
is( $active_parameter{array}[0],
    $expected_absolute_path, q{Updated array element to absolute path} );

## Then scalar path should now be absolute
is( $active_parameter{scalar},
    $expected_absolute_path, q{Updated scalar to absolute path} );

done_testing();
