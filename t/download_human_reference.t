#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp;
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2014 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_mip_hashes test_standard_cli };

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
        q{MIP::Recipes::Download::Human_reference} => [qw{ download_human_reference }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Download::Human_reference qw{ download_human_reference };

diag(   q{Test download_human_reference from Human_reference.pm v}
      . $MIP::Recipes::Download::Human_reference::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $test_dir  = File::Temp->newdir();
my $file_path = catfile( $test_dir, q{recipe_script.sh} );
my $log       = test_log( { log_name => q{MIP_DOWNLOAD}, } );

## Given
my %parameter = test_mip_hashes(
    {
        mip_hash_name => q{download_parameter},
    }
);
$parameter{project_id}    = q{test};
$parameter{reference_dir} = catfile($test_dir);
my %job_id;
my $recipe_name       = q{human_reference};
my $genome_version    = q{grch37};
my $reference_version = q{decoy_5};
my $reference_href    = $parameter{$recipe_name}{$genome_version}{$reference_version};
my $is_ok             = download_human_reference(
    {
        job_id_href    => \%job_id,
        parameter_href => \%parameter,
        reference_href => $reference_href,
        recipe_name    => $recipe_name,
        temp_directory => catfile($test_dir),
    }
);

## Then
ok( $is_ok, q{ Executed download recipe } . $recipe_name );

done_testing();
