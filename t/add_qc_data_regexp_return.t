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
use MIP::Constants qw{ $COMMA $SPACE };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Qc_data}        => [qw{ add_qc_data_regexp_return }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Qc_data qw{ add_qc_data_regexp_return };

diag(   q{Test add_qc_data_regexp_return from Qc_data.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a valid regexp when file exists
my $data_file_path = catfile( dirname($Bin), q{mip} );
my %qc_data;
my $recipe_name = q{mip};
my $regexp = q{perl -nae 'if ($_ =~ /\bour\s\$VERSION\b/xms) {print q{Got version};}'};
my $regexp_key = q{test_add_data_from_regexp};

add_qc_data_regexp_return(
    {
        data_file_path => $data_file_path,
        qc_href        => \%qc_data,
        recipe_name    => $recipe_name,
        regexp         => $regexp,
        regexp_key     => $regexp_key,
    }
);

my %expected_qc_data = ( $recipe_name => { $regexp_key => [ qw{Got version}, ], } );

## Then reg exp returned data should have been added to qc_data
is_deeply( \%qc_data, \%expected_qc_data, q{Added regexp return from system call} );

## Given a regexp when faulty return
my $faulty_regexp = q{perl -nae 'if ($_ =~ /\bnot our\s\$VERSION\b/xms) {}'};

my $return = add_qc_data_regexp_return(
    {
        data_file_path => $data_file_path,
        qc_href        => \%qc_data,
        recipe_name    => $recipe_name,
        regexp         => $faulty_regexp,
        regexp_key     => $regexp_key,
    }
);

## Then
is( $return, undef, q{Could not find a separator} );

done_testing();
