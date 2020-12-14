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
use MIP::Constants qw{ $COMMA $COLON $PIPE $SPACE };


BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Vcfparser}      => [qw{ write_feature_file_csq }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Vcfparser qw{ write_feature_file_csq };

diag(   q{Test write_feature_file_csq from Vcfparser.pm}
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

## Given feature files and transcripts
my $range_csq  = join $PIPE, qw{ a range_transcript };
my $select_csq = join $PIPE, qw{ a select_transcript };
my %vcf_record = (
    INFO_key_value     => { CSQ => 1, },
    range_transcripts  => [$range_csq],
    select_transcripts => [$select_csq],
);

my $info_field_counter = write_feature_file_csq(
    {
        filehandle         => $filehandle,
        info_field_counter => 0,
        select_fh          => $select_fh,
        vcf_record_href    => \%vcf_record,
    }
);

## Close the filehandle
close $filehandle;
close $select_fh;

## Then write feature file csq transcripts
my ($ret_range_csq) = $range_file_content =~ /\A CSQ=$range_csq \z/msx;
ok( $ret_range_csq, q{Wrote range csq transcripts} );
my ($ret_select_csq) = $select_file_content =~ /\A CSQ=$select_csq \z/msx;
ok( $ret_select_csq, q{Wrote select csq transcripts} );

## Then delete feature file CSQ record from input hash
is( $vcf_record{INFO_key_value}{CSQ},
    undef, q{Deleted feature file CSQ record from input hash} );

## Then increment info_field_counter
ok( $info_field_counter, q{Incremented info field counter} );

done_testing();
