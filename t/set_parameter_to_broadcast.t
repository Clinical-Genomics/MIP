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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.02;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $COMMA               => q{,};
Readonly my $SPACE               => q{ };
Readonly my $HASH_OF_HASH_INDEX  => q{3};
Readonly my $HASH_OF_ARRAY_INDEX => q{4};

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Set::Parameter} => [qw{ set_parameter_to_broadcast}],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Set::Parameter qw{ set_parameter_to_broadcast };

diag(   q{Test set_parameter_to_broadcast from Parameter.pm v}
      . $MIP::Set::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given different containers to broadcast
my %active_parameter = (
    Basil  => [qw{ Don't mention the war }],
    cast   => { starring => [qw{ John_Cleese Prunella_Scales Andrew_Sachs }], },
    Fawlty => { Towers => { genre => q{sitcom}, }, },
    Manuel => {
        line => q{Que},
    },
    Sybil => q{Ooohh, I knoooow},
);

# Add a order to the broadcast
my @order_parameters = qw{ Basil Manuel Sybil Fawlty cast };

# How to seperate elements in arrays
my %parameter = (
    Basil => {
        element_separator => $SPACE,
    }
);

## Store Broadcast message
my @broadcasts;

set_parameter_to_broadcast(
    {
        active_parameter_href => \%active_parameter,
        broadcasts_ref        => \@broadcasts,
        order_parameters_ref  => \@order_parameters,
        parameter_href        => \%parameter,
    }
);

## Then expect output according to container level
is( $broadcasts[0], q{Set Basil to: Don't mention the war}, q{Set array} );
is( $broadcasts[1], q{Set Manuel to: line=Que},             q{Set hash} );
is( $broadcasts[2], q{Set Sybil to: Ooohh, I knoooow},      q{Set scalar} );
is(
    $broadcasts[$HASH_OF_HASH_INDEX],
    q?Set Fawlty to: {Towers => genre=sitcom,} ?,
    q{Set hash of hash}
);

is(
    $broadcasts[$HASH_OF_ARRAY_INDEX],
    q?Set cast to: starring=John_Cleese Prunella_Scales Andrew_Sachs?,
    q{Set hash of array}
);

done_testing();
