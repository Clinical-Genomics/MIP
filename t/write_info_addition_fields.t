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
use MIP::Constants qw{ $COMMA $COLON $PIPE $SEMICOLON $SPACE };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Vcfparser}      => [qw{ write_info_addition_fields }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Vcfparser qw{ write_info_addition_fields };

diag(   q{Test write_info_addition_fields from Vcfparser.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

# For storing info to write
my $range_file_content;

## Store range feature file content in memory by using referenced variable
open my $filehandle, q{>}, \$range_file_content
  or croak q{Cannot write to}
  . $SPACE
  . $range_file_content
  . $COLON
  . $SPACE
  . $OS_ERROR;

# For storing info to write
my $select_file_content;

## Store select file content in memory by using referenced variable
open my $select_fh, q{>}, \$select_file_content
  or croak q{Cannot write to}
  . $SPACE
  . $select_file_content
  . $COLON
  . $SPACE
  . $OS_ERROR;

my $range_csq  = join $PIPE,      qw{ a range_transcript };
my $select_csq = join $PIPE,      qw{ a select_transcript };
my $info       = join $SEMICOLON, qw{ AF=0.4 IMPRECISE most_severe_pli=1 };
my %vcf_record = (
    INFO_addition => {
        AF              => 1,
        most_severe_pli => 1,
    },
    INFO_addition_select_feature => {
        SELECT_AF       => 1,
        most_severe_pli => 1,
    },
    INFO_addition_range_feature => {
        RANGE_AF        => 1,
        most_severe_pli => 1,
    },
    range_transcripts  => [$range_csq],
    select_transcripts => [$select_csq],
);

write_info_addition_fields(
    {
        filehandle      => $filehandle,
        select_fh       => $select_fh,
        vcf_record_href => \%vcf_record,
    }
);

## Close the filehandle
close $filehandle;
close $select_fh;

my $expected_info_add_content     = q{[;]AF=1};
my $expected_info_add_sel_content = q{[;]SELECT_AF=1};
my $expected_info_add_ran_content = q{[;]RANGE_AF=1};

## Then write info addition field to feature files
my ($ret_range_info) = $range_file_content =~ / $expected_info_add_content /msx;
ok( $ret_range_info, q{Wrote range info addition field to range file} );
my ($ret_select_info) = $select_file_content =~ / $expected_info_add_content /msx;
ok( $ret_select_info, q{Wrote select info addition field to select file} );

## Then write info addition select field to feature file
($ret_select_info) = $select_file_content =~ / $expected_info_add_sel_content /msx;
ok( $ret_select_info, q{Wrote select info addition select field to select file} );

## Then write info addition range field to feature file
($ret_range_info) = $range_file_content =~ / $expected_info_add_ran_content /msx;
ok( $ret_range_info, q{Wrote range info addition range field to range file} );

done_testing();
