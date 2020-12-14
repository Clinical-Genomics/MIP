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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Parameter}      => [qw{ get_parameter_attribute }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Parameter qw{ get_parameter_attribute };

diag(   q{Test get_parameter_attribute from Parameter.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a parameter when attribute is uninitilized
my $parameter_name = q{bwa_mem_file};
my %parameter      = (
    $parameter_name => {
        analysis_mode       => q{sample},
        associated_recipe   => [qw{mip bwa_mem}],
        build_file          => 1,
        data_type           => q{SCALAR},
        default             => { path_to_file => q{a_file}, },
        exists_check        => 1,
        file_tag            => q{build},
        outfile_suffix      => q{.sai},
        program_executables => [qw{ bwa samtools }],
        recipe_type         => q{aligners},
        type                => q{recipe},
    }
);

## When no attribute is supplied
my %attribute = get_parameter_attribute(
    {
        parameter_href => \%parameter,
        parameter_name => $parameter_name,
    }
);

## Then return entire hash
is_deeply( \%{ $parameter{$parameter_name} }, \%attribute, q{Got entire attribute hash} );

## When attribute is uninitilized
my $is_ok = get_parameter_attribute(
    {
        parameter_href => \%parameter,
        parameter_name => $parameter_name,
        attribute      => q{is_reference},
    }
);

## Then return undef
is( $is_ok, undef, q{Skipped attribute} );

## When scalar
my $is_scalar = get_parameter_attribute(
    {
        parameter_href => \%parameter,
        parameter_name => $parameter_name,
        attribute      => q{exists_check},
    }
);

## Then return scalar
ok( $is_scalar, q{Got scalar attribute} );

## When array
my @associated_recipes = get_parameter_attribute(
    {
        parameter_href => \%parameter,
        parameter_name => $parameter_name,
        attribute      => q{associated_recipe}
    }
);

## Then return array
is_deeply( \@{ $parameter{$parameter_name}{associated_recipe} },
    \@associated_recipes, q{Got array attribute} );

## When hash
my %parameter_default = get_parameter_attribute(
    {
        parameter_href => \%parameter,
        parameter_name => $parameter_name,
        attribute      => q{default},
    }
);

## Then return hash
is_deeply( \%{ $parameter{$parameter_name}{default} },
    \%parameter_default, q{Got hash attribute} );

done_testing();
