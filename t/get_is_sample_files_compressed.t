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
        q{MIP::File_info}      => [qw{ get_is_sample_files_compressed }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File_info qw{ get_is_sample_files_compressed };

diag(   q{Test get_is_sample_files_compressed from File_info.pm v}
      . $MIP::File_info::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a sample_id
my %file_info;
my $sample_id = q{a_sample_id};

## When no compression exists
my $return = get_is_sample_files_compressed(
    {
        file_info_href => \%file_info,
        sample_id      => $sample_id,
    }
);

## Then return undef
is( $return, undef, q{Return undef if not defined} );

## Given a true compression status
$file_info{is_files_compressed}{$sample_id} = 1;

## When all files are compressed
my $is_files_compressed = get_is_sample_files_compressed(
    {
        file_info_href => \%file_info,
        sample_id      => $sample_id,
    }
);

## Then return true
ok( $is_files_compressed, q{Return true if all files are compressed} );

## Given a false compression status
$file_info{is_files_compressed}{$sample_id} = 0;

## When all files are compressed
$is_files_compressed = get_is_sample_files_compressed(
    {
        file_info_href => \%file_info,
        sample_id      => $sample_id,
    }
);

## Then return false
is( $is_files_compressed, 0, q{Return false if not all files are compressed} );
done_testing();
