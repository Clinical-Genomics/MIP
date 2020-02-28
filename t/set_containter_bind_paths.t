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
use MIP::Constants qw{ $BACKWARD_SLASH $COLON $COMMA $DOUBLE_QUOTE $SPACE };
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
        q{MIP::Set::Parameter} => [qw{ set_container_bind_paths }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Set::Parameter qw{ set_container_bind_paths };

diag(   q{Test set_container_bind_paths from Parameter.pm v}
      . $MIP::Set::Parameter::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a container without bind paths
my %container = (
    executable => {
        binary => undef,
    },
    uri => q{docker address},
);
my @bind_paths = ( catdir(qw{ a path }) );
set_container_bind_paths(
    {
        bind_paths_ref => \@bind_paths,
        container_href => \%container,
    }
);

## Then create key and set path
is_deeply( $container{program_bind_paths},
    \@bind_paths, q{Create and set program_bind_paths} );

## Given an already populated program_bind_paths key
my @extra_paths = ( catdir(qw{ another path }) );
set_container_bind_paths(
    {
        bind_paths_ref => \@extra_paths,
        container_href => \%container,
    }
);

## Then add to array
my @expected = ( catdir(qw{ a path }), catdir(qw{ another path }) );
is_deeply( $container{program_bind_paths}, \@expected, q{Append to program_bind_paths} );

done_testing();
