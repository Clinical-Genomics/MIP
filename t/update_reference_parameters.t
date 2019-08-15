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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Test::Fixtures qw{ test_standard_cli };

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
        q{MIP::Update::Parameters} => [qw{ update_reference_parameters }],
        q{MIP::Test::Fixtures}     => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Update::Parameters qw{ update_reference_parameters };

diag(   q{Test update_reference_parameters from Parameters.pm v}
      . $MIP::Update::Parameters::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $reference_dir = catdir(qw{ test ref_dir });

my %parameter = (
    array       => { associated_recipe => [qw{ recipe_1 }], },
    file_name_1 => { associated_recipe => [qw{ recipe_1 }], },
    file_name_2 => { associated_recipe => [qw{ recipe_2 }], },
    hash        => { associated_recipe => [qw{ recipe_1 }], },
);

my %active_parameter = (
    array         => [qw{ file_1 file_2 }],
    file_name_1   => q{file_0},
    file_name_2   => q{file_2},
    hash          => { file_3 => q{info_key_1} },
    recipe_1      => 1,
    reference_dir => $reference_dir,
);

PARAMETER:
foreach my $parameter_name ( keys %parameter ) {

    update_reference_parameters(
        {
            active_parameter_href => \%active_parameter,
            associated_recipes_ref =>
              \@{ $parameter{$parameter_name}{associated_recipe} },
            parameter_name => $parameter_name,
        }
    );
}

is( $active_parameter{file_name_2}, q{file_2}, q{Skipped setting file reference path} );

is(
    $active_parameter{file_name_1},
    catfile( $reference_dir, q{file_0} ),
    q{Set file reference path}
);

is(
    $active_parameter{array}[0],
    catfile( $reference_dir, q{file_1} ),
    q{Set array reference path}
);

UPDATED_FILE:
foreach my $updated_file ( keys %{ $active_parameter{hash} } ) {

    is( $updated_file, catfile( $reference_dir, q{file_3} ), q{Set hash reference path} );

}

done_testing();
