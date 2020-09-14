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

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $DOT $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.02;

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
        q{MIP::File::Path}     => [qw{ remove_file_path_suffix }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File::Path qw{ remove_file_path_suffix };

diag(   q{Test remove_file_path_suffix from Path.pm v}
      . $MIP::File::Path::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a file name path with a suffix
my $file_path = catfile(qw{ a dir file.fastq });

## When supplied suffix does not match target file name suffix
my $outfile_path_no_suffix = remove_file_path_suffix(
    {
        file_path         => $file_path,
        file_suffixes_ref => [ $DOT . q{gz} ],
    }
);

## Then return undef
is( $outfile_path_no_suffix, undef, q{Did not match suffix(es)} );

## Given a file name path with a suffix
my $matching_suffix_file_path = catfile(qw{ a dir file.fastq.gz });

## When suffix match file name
$outfile_path_no_suffix = remove_file_path_suffix(
    {
        file_path         => $matching_suffix_file_path,
        file_suffixes_ref => [ $DOT . q{gz} ],
    }
);

## Then return file path without suffix
is( $outfile_path_no_suffix, $file_path, q{Match suffix} );

## Given a file name path with 2 suffixes
my $matching_long_suffix_file_path = catfile(qw{ a dir files.tar.gz });

## When double suffix match target file name suffix
$outfile_path_no_suffix = remove_file_path_suffix(
    {
        file_path         => $matching_long_suffix_file_path,
        file_suffixes_ref => [ $DOT . q{tar} . $DOT . q{gz} ],
    }
);

$file_path = catfile(qw{ a dir files });

## Then return file path without suffixes
is( $outfile_path_no_suffix, $file_path, q{Match partial suffix} );

## Given a file name
my $compressed_file_name = catfile(qw{ files.tar.gz });

## When no dir and matching suffix
$outfile_path_no_suffix = remove_file_path_suffix(
    {
        file_path         => $compressed_file_name,
        file_suffixes_ref => [ $DOT . q{gz} ],
    }
);

my $file_name = q{files.tar};

## Then return file name without suffix
is( $outfile_path_no_suffix, $file_name, q{Match suffix no dir} );

done_testing();
