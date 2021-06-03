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
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Active_parameter} => [qw{ set_vcfparser_outfile_counter }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Active_parameter qw{ set_vcfparser_outfile_counter };

diag(   q{Test set_vcfparser_outfile_counter from Active_parameter.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a vcfparser select file
my %active_parameter_test = (
    sv_vcfparser          => { type => q{recipe} },
    vcfparser_ar          => { type => q{recipe} },
    vcfparser_select_file => 1,
);

## When vcfparser is used with select file and sv_vcfparser without.
set_vcfparser_outfile_counter( { active_parameter_href => \%active_parameter_test, } );

## Then vcfparser_outfile_count should be 2
is( $active_parameter_test{vcfparser_outfile_count},
    2, q{vcfparser_ar used with a select file -> 2 outfiles} );

## Then sv_vcfparser_outfile_count should be 1
is( $active_parameter_test{sv_vcfparser_outfile_count},
    1, q{sv_vcfparser used without a select file -> 1 outfile} );

## Given a vcfparser select file and sv vcfparser select file
%active_parameter_test = (
    sv_vcfparser             => { type => q{recipe} },
    vcfparser_ar             => { type => q{recipe} },
    sv_vcfparser_select_file => 1,
    vcfparser_select_file    => 1,
);

## When both vcfparser and sv_vcfparser are used with select files.
set_vcfparser_outfile_counter( { active_parameter_href => \%active_parameter_test, } );

## Then vcfparser_outfile_count should be 2
is( $active_parameter_test{vcfparser_outfile_count},
    2, q{vcfparser_ar used with a select file -> 2 outfiles} );

## Then sv_vcfparser_outfile_count should be 2
is( $active_parameter_test{sv_vcfparser_outfile_count},
    2, q{sv_vcfparser used with a select file -> 2 outfiles} );

done_testing();
