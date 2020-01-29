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
use MIP::Constants qw{ $COMMA $SPACE };
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
        q{MIP::Sample_info}    => [qw{ set_file_path_to_store }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Sample_info qw{ set_file_path_to_store };

diag(   q{Test set_file_path_to_store from Sample_info.pm v}
      . $MIP::Sample_info::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given file info
my $format      = q{vcf};
my $path        = catfile(qw{ a test path.vcf.gz });
my $path_index  = catfile(qw{ a test path.vcf.gz.tbi });
my $recipe_name = q{gatk_variantrecalibration};
my $sample_id   = q{sample_1};
my %sample_info;
my $tag = q{a_tag};

set_file_path_to_store(
    {
        format           => $format,
        id               => $sample_id,
        path             => $path,
        path_index       => $path_index,
        recipe_name      => $recipe_name,
        sample_info_href => \%sample_info,
        tag              => $tag,
    }
);

my %expected_sample_info = (
    files => [
        {
            format     => $format,
            id         => $sample_id,
            path       => $path,
            path_index => $path_index,
            step       => $recipe_name,
            tag        => $tag,
        },
    ],
);

## Then file path should be set
is_deeply( \%sample_info, \%expected_sample_info, q{Set file path in hash} );

done_testing();
