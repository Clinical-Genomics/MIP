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
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $COMMA      => q{,};
Readonly my $DOT        => q{.};
Readonly my $SPACE      => q{ };
Readonly my $UNDERSCORE => q{_};

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Sample_info}    => [qw{ set_most_complete_vcf }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Sample_info qw{ set_most_complete_vcf };

diag(   q{Test set_most_complete_vcf from Sample_info.pm v}
      . $MIP::Sample_info::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $path             = catfile(qw{ a test vcf_file });
my $recipe_name_test = q{sv_rankvariants};
my $file_suffix      = $DOT . q{vcf};
my $vcf_file_key =
  q{sv} . $UNDERSCORE . substr( $file_suffix, 1 ) . $UNDERSCORE . q{file};

my %active_parameter = ( $recipe_name_test => 1, );

my %sample_info;
set_most_complete_vcf(
    {
        active_parameter_href     => \%active_parameter,
        path                      => $path,
        recipe_name               => $recipe_name_test,
        sample_info_href          => \%sample_info,
        vcfparser_outfile_counter => 0,
        vcf_file_key              => $vcf_file_key,
    }
);
is( $sample_info{$vcf_file_key}{research}{path}, $path, q{Added research path} );

set_most_complete_vcf(
    {
        active_parameter_href     => \%active_parameter,
        path                      => $path,
        recipe_name               => $recipe_name_test,
        sample_info_href          => \%sample_info,
        vcfparser_outfile_counter => 1,
        vcf_file_key              => $vcf_file_key,
    }
);
is( $sample_info{$vcf_file_key}{clinical}{path}, $path, q{Added clinical path} );

done_testing();
