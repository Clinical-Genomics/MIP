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
        q{MIP::File_info}      => [qw{ add_sample_no_direction_infile_prefixes }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File_info qw{ add_sample_no_direction_infile_prefixes };

diag(   q{Test add_sample_no_direction_infile_prefixes from File_info.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a sample_id and an no_direction_infile_prefixes
my %file_info;
my $no_direction_infile_prefixes = q{ADM1059A1_161011_HHJJCCCXY_NAATGCGC_lane7};
my $sample_id                  = q{sample_id};

## Then do not set no_direction_infile_prefixes in file_info hash
add_sample_no_direction_infile_prefixes(
    {
        file_info_href  => \%file_info,
        mip_file_format => $no_direction_infile_prefixes,
        sample_id       => $sample_id,
    }
);

## Then set no_direction_infile_prefixes in file_info hash
is( @{ $file_info{$sample_id}{no_direction_infile_prefixes} },
    1, q{Add no_direction_infile_prefixes to file_info} );

done_testing();
