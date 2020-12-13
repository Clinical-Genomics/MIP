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
        q{MIP::File_info}      => [qw{ set_file_tag }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File_info qw{ set_file_tag };

diag(   q{Test set_file_tag from File_info.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given id, recipe and file_tag
my %file_info;
my $file_tag    = q{_sort};
my $id          = q{a_sample_id};
my $recipe_name = q{bwa_mem};

set_file_tag(
    {
        file_info_href => \%file_info,
        file_tag       => $file_tag,
        id             => $id,
        recipe_name    => $recipe_name,
    }
);

## Then file tag should be set for id and recipe
is( $file_info{$id}{$recipe_name}{file_tag},
    $file_tag, q{Set file tag for id and recipe} );

done_testing();
