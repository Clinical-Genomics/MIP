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
        q{MIP::Active_parameter} => [qw{ set_exome_target_bed }],

    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Active_parameter qw{ set_exome_target_bed };

diag(   q{Test set_exome_target_bed from Active_parameter.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a bed file and a sample_id string
my %active_parameter      = ();
my $exome_target_bed_file = q{a_bed_file};
my $sample_id_string      = join $COMMA, qw{ sample_1 sample_2};

set_exome_target_bed(
    {
        active_parameter_href => \%active_parameter,
        exome_target_bed_file => $exome_target_bed_file,
        sample_id_string      => $sample_id_string,
    }
);

## Then sample_ids string should be attached to exome target bed
is( $active_parameter{exome_target_bed}{$exome_target_bed_file},
    q{sample_1,sample_2}, q{Added sample_ids string to exome target bed} );

done_testing();
