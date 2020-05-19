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
        q{MIP::File_info}      => [qw{ add_sample_infile_both_strands_prefix }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File_info qw{ add_sample_infile_both_strands_prefix };

diag(   q{Test add_sample_infile_both_strands_prefix from File_info.pm v}
      . $MIP::File_info::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a sample_id
my %file_info;
my $mip_file_format_with_direction = q{a_file_format_1};
my $sample_id                      = q{sample_id};

add_sample_infile_both_strands_prefix(
    {
        file_info_href                 => \%file_info,
        mip_file_format_with_direction => $mip_file_format_with_direction,
        sample_id                      => $sample_id,
    }
);

## Then add infile_both_strands_prefix
is( @{ $file_info{$sample_id}{infile_both_strands_prefix} },
    1, q{Added infile_both_strands_prefix to file_info} );

done_testing();
