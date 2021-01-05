#! /usr/bin/env perl

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
use MIP::Constants qw{ $COMMA $SPACE $UNDERSCORE };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = ( q{MIP::File_info} => [qw{ set_merged_infile_prefix }], );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File_info qw{ set_merged_infile_prefix };

diag(   q{Test set_merged_infile_prefix from File_info.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a lanes id and sample id
my %file_info;
my $lanes_id             = q{123};
my $sample_id            = q{sample_1};
my $merged_infile_prefix = $sample_id . $UNDERSCORE . q{lanes} . $UNDERSCORE . $lanes_id;

## When setting the merge_infile_prefix
set_merged_infile_prefix(
    {
        file_info_href       => \%file_info,
        merged_infile_prefix => $merged_infile_prefix,
        sample_id            => $sample_id,
    }
);

my $set_merged_infile_prefix = $file_info{$sample_id}{merged_infile};

## Then the merge infile prefix should be set in the file_info hash
is( $set_merged_infile_prefix, $merged_infile_prefix, q{Set merged infile prefix for sample id} );

done_testing();
