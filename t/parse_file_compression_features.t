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
        q{MIP::File_info}      => [qw{ parse_file_compression_features }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File_info qw{ parse_file_compression_features };

diag(   q{Test parse_file_compression_features from File_info.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given uncompressed file
my %file_info;
my $file_name = q{a_file.fastq};
my $sample_id = q{a_sample_id};

my $read_file_command = parse_file_compression_features(
    {
        file_info_href => \%file_info,
        file_name      => $file_name,
        sample_id      => $sample_id,
    }
);

## Then return false
is( $file_info{$sample_id}{$file_name}{is_file_compressed},
    0, q{File was not compressed} );

## Then set read command to handle uncompressed file
is( $file_info{$sample_id}{$file_name}{read_file_command},
    q{cat}, q{Set read file command for uncompressed file} );

## Then returned read command should be for uncompressed files
is( $read_file_command, q{cat}, q{Returned read file command for uncompressed files} );

## Given compressed file
$file_name .= q{a_file.fastq.gz};

$read_file_command = parse_file_compression_features(
    {
        file_info_href => \%file_info,
        file_name      => $file_name,
        sample_id      => $sample_id,
    }
);

## Then return true
is( $file_info{$sample_id}{$file_name}{is_file_compressed}, 1, q{File was compressed} );

## Then set read command to handle uncompressed file
is( $file_info{$sample_id}{$file_name}{read_file_command},
    q{gzip -d -c}, q{Set read file command for compressed file} );

## Then returned read command should be for uncompressed files
is( $read_file_command, q{gzip -d -c},
    q{Returned read file command for compressed files} );

done_testing();
