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
use Modern::Perl qw{ 2018 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COLON $COMMA $SPACE };
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
        q{MIP::Recipes::Download::Get_reference} => [qw{ get_reference }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Recipes::Download::Get_reference qw{ get_reference };

diag(   q{Test get_reference from Get_reference.pm v}
      . $MIP::Recipes::Download::Get_reference::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

# Create anonymous filehandle
my $FILEHANDLE = IO::Handle->new();

# For storing info to write
my $file_content;

## Store file content in memory by using referenced variable
open $FILEHANDLE, q{>}, \$file_content
  or croak q{Cannot write to} . $SPACE . $file_content . $COLON . $SPACE . $OS_ERROR;

my $test_dir = File::Temp->newdir();
my $log      = test_log( { log_name => q{MIP_DOWNLOAD}, no_screen => 1, } );

## Given analysis parameters
my $genome_version    = q{grch37};
my $recipe_name       = q{human_reference};
my $reference_version = q{decoy_5};

my %active_parameter = test_mip_hashes(
    {
        mip_hash_name => q{download_active_parameter},
    }
);
$active_parameter{project_id}    = q{test};
$active_parameter{reference_dir} = catfile($test_dir);
my $reference_href =
  $active_parameter{reference_feature}{$recipe_name}{$genome_version}{$reference_version};

my $is_ok = get_reference(
    {
        FILEHANDLE     => $FILEHANDLE,
        recipe_name    => $recipe_name,
        reference_dir  => $active_parameter{reference_dir},
        reference_href => $reference_href,
    }
);

## Then
ok( $is_ok, q{ Executed download recipe } . $recipe_name );

close $FILEHANDLE;
done_testing();
