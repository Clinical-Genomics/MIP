#! /usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile  };
use File::Spec::Functions qw{ catdir catfile  };
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
        q{MIP::Language::Perl} => [qw{ get_cpan_file_modules }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Language::Perl qw{ get_cpan_file_modules };

diag(   q{Test get_cpan_file_modules from Perl.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a cpanfile
my $cpan_file_path = catfile( $Bin, qw{ data test_data cpanfile } );

## When reading cpanfile to get modules
my @cpan_modules = get_cpan_file_modules( { cpanfile_path => $cpan_file_path, } );

my @expected_modules = qw{ Array::Utils Clone };

## Then all the cpan modules in cpan file path should be returned in lexiographical order
is_deeply( \@cpan_modules, \@expected_modules, q{ Got cpan module(s)} );

done_testing();
