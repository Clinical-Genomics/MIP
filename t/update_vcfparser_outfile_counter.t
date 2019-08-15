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
Readonly my $COMMA => q{,};
Readonly my $SPACE => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Update::Parameters} => [qw{ update_vcfparser_outfile_counter }],
        q{MIP::Test::Fixtures}     => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Update::Parameters qw{ update_vcfparser_outfile_counter };

diag(   q{Test update_vcfparser_outfile_counter from Parameters.pm v}
      . $MIP::Update::Parameters::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

# Test the number of oufiles when vcfparser is used with select file and sv_vcfparser without.
my %active_parameter_test = (
    sv_vcfparser          => { type => q{recipe} },
    vcfparser_ar          => { type => q{recipe} },
    vcfparser_select_file => 1,
);

update_vcfparser_outfile_counter( { active_parameter_href => \%active_parameter_test, } );

is( $active_parameter_test{vcfparser_outfile_count},
    2, q{vcfparser_ar used with a select file -> 2 outfiles. Test passed.} );
is( $active_parameter_test{sv_vcfparser_outfile_count},
    1, q{sv_vcfparser used without a select file -> 1 outfile. Test passed.} );

# Test the number of oufiles when both vcfparser and sv_vcfparser are used with select files.
%active_parameter_test = (
    sv_vcfparser             => { type => q{recipe} },
    vcfparser_ar             => { type => q{recipe} },
    sv_vcfparser_select_file => 1,
    vcfparser_select_file    => 1,
);

update_vcfparser_outfile_counter( { active_parameter_href => \%active_parameter_test, } );

is( $active_parameter_test{vcfparser_outfile_count},
    2, q{vcfparser_ar used with a select file -> 2 outfiles. Test passed.} );
is( $active_parameter_test{sv_vcfparser_outfile_count},
    2, q{sv_vcfparser used with a select file -> 2 outfiles. Test passed.} );

done_testing();
