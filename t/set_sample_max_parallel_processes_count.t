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
use MIP::Test::Fixtures qw{ test_standard_cli };

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
        q{MIP::File_info}      => [qw{ set_sample_max_parallel_processes_count }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File_info qw{ set_sample_max_parallel_processes_count };

diag(   q{Test set_sample_max_parallel_processes_count from File_info.pm v}
      . $MIP::File_info::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a sample id
my %file_info;
my $sample_id = q{ADM1059A1};

## When max_parallel processing counter is 0
my $max_parallel_processes_count = 0;

set_sample_max_parallel_processes_count(
    {
        file_info_href               => \%file_info,
        max_parallel_processes_count => $max_parallel_processes_count,
        sample_id                    => $sample_id,
    }
);

## Then set max_parallel_processes_count to 0
is( $file_info{max_parallel_processes_count}{$sample_id},
    0, q{Set max_parallel_processes_count per sample id} );

done_testing();
