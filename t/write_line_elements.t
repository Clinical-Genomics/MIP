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
use MIP::Constants qw{ $COMMA $COLON $DOT $NEWLINE $PIPE $SPACE $TAB };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Vcfparser}      => [qw{ write_line_elements }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Vcfparser qw{ write_line_elements };

diag(   q{Test write_line_elements from Vcfparser.pm}
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

my $first_separator = $TAB;
my $last_index      = 2;
my $last_separator  = $NEWLINE;
my @line_elements   = ( 1, 2, $DOT, q{PASS} );
my $start_index     = 0;
my $range_csq       = join $PIPE, qw{ a range_transcript };
my $select_csq      = join $PIPE, qw{ a select_transcript };
my %vcf_record      = (
    range_transcripts  => [$range_csq],
    select_transcripts => [$select_csq],
);

write_line_elements(
    {
        filehandle        => $filehandle,
        first_separator   => $first_separator,
        last_index        => $last_index,
        last_separator    => $last_separator,
        line_elements_ref => \@line_elements,
        select_fh         => $select_fh,
        start_index       => 0,
        vcf_record_href   => \%vcf_record,
    }
);

## Close the filehandle
close $filehandle;
close $select_fh;

my $expected_line = q{\s+ 1 \s+ 2 \s+ } . $DOT . $NEWLINE;

## Then write line elements for feature files
my ($ret_range_info) = $range_file_content =~ /\A $expected_line  /sxm;
ok( $ret_range_info, q{Wrote range line elements} );
my ($ret_select_info) = $select_file_content =~ /\A $expected_line /msx;
ok( $ret_select_info, q{Wrote select line elements} );

done_testing();
