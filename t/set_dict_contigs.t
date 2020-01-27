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
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_log test_standard_cli };

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
        q{MIP::File_info}      => [qw{ set_dict_contigs }],
        q{MIP::Parameter}      => [qw{ set_parameter_build_file_status }],
        q{MIP::Test::Fixtures} => [qw{ test_log test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File_info qw{ set_dict_contigs };
use MIP::Parameter qw{ set_parameter_build_file_status };

diag(   q{Test set_dict_contigs from Reference.pm v}
      . $MIP::File_info::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Creates log object
my $log = test_log( {} );

## Given proper data
my %file_info;
my $dict_file_path = catfile( $Bin, qw{ data references grch37_homo_sapiens_-d5-.dict } );
my %parameter = ( human_genome_reference_file_endings => { build_file => 0, } );

set_dict_contigs(
    {
        dict_file_path => $dict_file_path,
        file_info_href => \%file_info,
        parameter_href => \%parameter,
    }
);

my @expected_contigs =
  qw{ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT GL000207.1 GL000226.1 GL000229.1 GL000231.1 GL000210.1 GL000239.1 GL000235.1 GL000201.1 GL000247.1 GL000245.1 GL000197.1 GL000203.1 GL000246.1 GL000249.1 GL000196.1 GL000248.1 GL000244.1 GL000238.1 GL000202.1 GL000234.1 GL000232.1 GL000206.1 GL000240.1 GL000236.1 GL000241.1 GL000243.1 GL000242.1 GL000230.1 GL000237.1 GL000233.1 GL000204.1 GL000198.1 GL000208.1 GL000191.1 GL000227.1 GL000228.1 GL000214.1 GL000221.1 GL000209.1 GL000218.1 GL000220.1 GL000213.1 GL000211.1 GL000199.1 GL000217.1 GL000216.1 GL000215.1 GL000205.1 GL000219.1 GL000224.1 GL000223.1 GL000195.1 GL000212.1 GL000222.1 GL000200.1 GL000193.1 GL000194.1 GL000225.1 GL000192.1 NC_007605 hs37d5 };

## Then return the dict contigs
is_deeply( \@expected_contigs, \@{ $file_info{contigs} }, q{Set dict file contigs} );

## Given a file when it needs to be built
set_parameter_build_file_status {
    (
        parameter_href => \%parameter,
        parameter_name => q{human_genome_reference_file_endings},
        status         => 1,
    )
};

my $is_ok = set_dict_contigs(
    {
        dict_file_path => $dict_file_path,
        file_info_href => \%file_info,
        parameter_href => \%parameter,
    }
);

## Then exit and throw FATAL log message
is( $is_ok, undef, q{Exit if file need to be built} );

done_testing();
