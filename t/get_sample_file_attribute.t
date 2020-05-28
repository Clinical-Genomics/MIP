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
use Test::Trap;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.01;

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
        q{MIP::File_info}      => [qw{ get_sample_file_attribute }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File_info qw{ get_sample_file_attribute };

diag(   q{Test get_sample_file_attribute from File_info.pm v}
      . $MIP::File_info::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given a sample with info
my $file_name     = q{a_file.gz};
my $infile_prefix = q{a_file};
my $sample_id     = q{a_sample_id};
my %file_info     = (
    $sample_id => {
        $file_name => {
            is_file_compressed => 1,
            read_file_command  => q{gzip -d -c},
        },
        no_direction_infile_prefixes => [ $infile_prefix, ],
    },
);
my %attribute = (
    is_file_compressed => 1,
    read_file_command  => q{gzip -d -c},
);

while ( my ( $attribute, $attribute_value ) = each %attribute ) {

    my $got_attribute = get_sample_file_attribute(
        {
            attribute      => $attribute,
            file_info_href => \%file_info,
            file_name      => $file_name,
            sample_id      => $sample_id,
        }
    );

    ## Then we should get an attribute
    is( $got_attribute, $attribute_value,
            qq{Got $sample_id file name: $file_name attribute: }
          . $attribute . q{ => }
          . $attribute_value );
}

## Given a no attribute and file name in call
my %got_sample_href = get_sample_file_attribute(
    {
        file_info_href => \%file_info,
        sample_id      => $sample_id,
    }
);

## Then return entire file info sample id hash
is_deeply(
    \%got_sample_href,
    \%{ $file_info{$sample_id} },
    q{Returned file info for sample id }
);

## Given a no attribute in call
## When file name point to hash
my %got_attribute_href = get_sample_file_attribute(
    {
        file_info_href => \%file_info,
        file_name      => $file_name,
        sample_id      => $sample_id,
    }
);

## Then return entire sample id file name attribute hash
is_deeply( \%got_attribute_href, \%attribute,
    q{Returned sample id file name attribute hash} );

## When file name points to array
my @got_attributes_ref = get_sample_file_attribute(
    {
        file_info_href => \%file_info,
        file_name      => q{no_direction_infile_prefixes},
        sample_id      => $sample_id,
    }
);

my @no_direction_infile_prefixes_attributes = ($infile_prefix);

## Then return entire sample id file name attribute attribute
is_deeply(
    \@got_attributes_ref,
    \@no_direction_infile_prefixes_attributes,
    q{Returned sample id file name attribute array}
);

## Given a undefined attribute in file_info hash
delete $file_info{$sample_id}{$file_name}{is_file_compressed};

my $got_attribute = get_sample_file_attribute(
    {
        attribute      => q{is_file_compressed},
        file_info_href => \%file_info,
        file_name      => $file_name,
        sample_id      => $sample_id,
    }
);

## Then return false
is( $got_attribute, undef, q{Returned undef for undefined sample file name attribute} );

done_testing();
