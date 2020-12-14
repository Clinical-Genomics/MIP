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
        q{MIP::Sample_info}    => [qw{ set_read_direction_file_attributes }],
);

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Sample_info qw{ set_read_direction_file_attributes };

diag(   q{Test set_read_direction_file_attributes from Sample_info.pm}
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a sample_id, file names and file attributes
my %both_directions_metric         = ( flowcell => q{a_flowcell}, );
my $mip_file_format_with_direction = q{a_file_name_1.fastq.gz};
my $mip_file_format                = q{a_file_name};
my $sample_id                      = q{a_sample_id};
my %sample_info;

set_read_direction_file_attributes(
    {
        direction_file_name => $mip_file_format_with_direction,
        file_name           => $mip_file_format,
        read_direction_href => \%both_directions_metric,
        sample_id           => $sample_id,
        sample_info_href    => \%sample_info,
    }
);

## Unpack
my $sample_info_file_attribute =
  $sample_info{sample}{$sample_id}{file}{$mip_file_format}{read_direction_file}
  {$mip_file_format_with_direction}{flowcell};

my $expected_sample_info_file_attribute = q{a_flowcell};

## Then
is(
    $sample_info_file_attribute,
    $expected_sample_info_file_attribute,
    q{Set read direction file attribute}
);

done_testing();
