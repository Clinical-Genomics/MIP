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
        q{MIP::File_info}      => [qw{ parse_files_compression_status }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File_info qw{ parse_files_compression_status };

diag(   q{Test parse_files_compression_status from File_info.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given uncompressed file
my $not_compressed_file = q{a_file.fastq};
my $sample_id           = q{a_sample_id};
my %file_info           = (
    $sample_id => {
        mip_infiles => [ $not_compressed_file, ],
        $not_compressed_file => { is_file_compressed => 0, }
    },
);

parse_files_compression_status(
    {
        file_info_href => \%file_info,
        sample_id      => $sample_id,
    }
);

## Then return false
is( $file_info{is_files_compressed}{$sample_id}, 0, q{All files were not compressed} );

## Given compressed file
my $compressed_file = q{a_file.fastq.gz};
push @{ $file_info{$sample_id}{mip_infiles} }, $compressed_file;
$file_info{$sample_id}{$compressed_file}{is_file_compressed} = 1;

parse_files_compression_status(
    {
        file_info_href => \%file_info,
        sample_id      => $sample_id,
    }
);

## Then return false
is( $file_info{is_files_compressed}{$sample_id}, 0, q{Not all files were compressed} );

## Given only a compressed file
$file_info{$sample_id}{mip_infiles} = [$compressed_file];
delete $file_info{$sample_id}{$not_compressed_file};

parse_files_compression_status(
    {
        file_info_href => \%file_info,
        sample_id      => $sample_id,
    }
);

## Then return true
ok( $file_info{is_files_compressed}{$sample_id}, q{All files were compressed} );
done_testing();
