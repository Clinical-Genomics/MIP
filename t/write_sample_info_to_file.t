#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Path qw{ remove_tree };
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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Sample_info}    => [qw{ write_sample_info_to_file }],
        q{MIP::Test::Fixtures} => [qw{ test_log }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Sample_info qw{ write_sample_info_to_file };

diag(   q{Test write_sample_info_to_file from Sample_info.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $log = test_log( { no_screen => 1, } );

## Given no sample info file to write to
my $sample_info_file;
my %sample_info;

my $is_ok = write_sample_info_to_file(
    {
        sample_info_file => $sample_info_file,
        sample_info_href => \%sample_info,
    }
);

## Then return false
is( $is_ok, undef, q{Do not write file} );

## Given a sample info file to write to
$sample_info_file = catfile( $Bin, q{a_test_write_sample_info_file.yaml} );

write_sample_info_to_file(
    {
        sample_info_file => $sample_info_file,
        sample_info_href => \%sample_info,
    }
);

ok( -e $sample_info_file, q{Wrote sample info file to disk} );

## Clean-up
remove_tree($sample_info_file);

done_testing();
