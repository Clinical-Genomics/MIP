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
use MIP::Constants qw{ $COMMA $DOT $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.01;

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
        q{MIP::Parse::File}    => [qw{ parse_file_suffix }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Parse::File qw{ parse_file_suffix };

diag(   q{Test parse_file_suffix from File.pm v}
      . $MIP::Parse::File::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a decompressed file when suffix does not match file name
my $file_name = catfile(qw{ a dir file.fastq });

my $outfile_path_no_suffix = parse_file_suffix(
    {
        file_name   => $file_name,
        file_suffix => $DOT . q{gz},
    }
);

## Then return undef
is( $outfile_path_no_suffix, undef, q{Did not match suffix} );

## Given a compressed file when suffix match file name
my $compressed_file_name = catfile(qw{ a dir file.fastq.gz });

$outfile_path_no_suffix = parse_file_suffix(
    {
        file_name   => $compressed_file_name,
        file_suffix => $DOT . q{gz},
    }
);

## Then return file name without suffix
is( $outfile_path_no_suffix, $file_name, q{Match suffix} );

## Given a compressed file when suffix match file name
$compressed_file_name = catfile(qw{ a dir files.tar.gz });

$outfile_path_no_suffix = parse_file_suffix(
    {
        file_name   => $compressed_file_name,
        file_suffix => $DOT . q{gz},
    }
);

$file_name = catfile(qw{ a dir files.tar });
## Then return file name without suffix
is( $outfile_path_no_suffix, $file_name, q{Match partial suffix} );
done_testing();
