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
use MIP::Constants qw{ $COMMA $SPACE };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Test::Commands} => [qw{ build_call }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Test::Commands qw{ build_call };

diag(   q{Test build_call from Commands.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given scalar input (input_value)
my $argument    = q{jedi};
my $input_value = q{luke};
my %required_argument =
  ( darth_vader => { input_value_href => { is_a => q{sith_lord}, }, }, );

my @scalar_args = build_call(
    {
        argument               => $argument,
        input_value            => $input_value,
        required_argument_href => \%required_argument,
    }
);
my @expected_scalar_args = qw{ darth_vader sith_lord jedi luke };

## Then built call should be returned
is_deeply( \@scalar_args, \@expected_scalar_args, q{Built scalar args} );

## Given array input (input_values_ref)
my @input_values = qw{ luke obi-wan };

my @array_args = build_call(
    {
        argument               => $argument,
        input_values_ref       => \@input_values,
        required_argument_href => \%required_argument,
    }
);

my @expected_array_args = ( qw{ darth_vader sith_lord jedi }, \@input_values );

## Then built call should be returned with array ref
is_deeply( \@array_args, \@expected_array_args, q{Built array args} );

## Given hash input (input_value_href)
my %input_value_hash = (
    luke       => q{padawan},
    q{obi-wan} => q{master},
);

my @hash_args = build_call(
    {
        argument               => $argument,
        input_value_href       => \%input_value_hash,
        required_argument_href => \%required_argument,
    }
);

my @expected_hash_args = ( qw{ darth_vader sith_lord jedi }, \%input_value_hash );

## Then built call should be returned with hash_ref
is_deeply( \@hash_args, \@expected_hash_args, q{Built hash args} );

## Given required hash input (input_value_href)
%input_value_hash = ( darth_vader => { is_a => q{sith_lord}, }, );

@hash_args = build_call(
    {
        argument               => q{darth_vader},
        input_value_href       => \%input_value_hash,
        required_argument_href => \%required_argument,
    }
);

@expected_hash_args = ( qw{ darth_vader sith_lord darth_vader }, \%input_value_hash );

## Then built call should be returned with hash ref
is_deeply( \@hash_args, \@expected_hash_args, q{Built hash args for required hash} );

done_testing();
