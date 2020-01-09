#!/usr/bin/env perl

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp;
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
use MIP::Test::Fixtures qw{ test_mip_hashes test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.00;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
## Constants
Readonly my $COMMA        => q{,};
Readonly my $COLON        => q{:};
Readonly my $DOT          => q{.};
Readonly my $DOUBLE_QUOTE => q{"};
Readonly my $NEWLINE      => qq{\n};
Readonly my $SINGLE_QUOTE => q{'};
Readonly my $SPACE        => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::File::Format::Yaml} => [qw{ write_yaml }],
        q{MIP::Test::Fixtures}     => [qw{ test_mip_hashes test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File::Format::Yaml qw{ write_yaml };

diag(   q{Test write_yaml from Yaml.pm v}
      . $MIP::File::Format::Yaml::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $test_dir = File::Temp->newdir();

## Given a hash
my %active_parameter = test_mip_hashes( { mip_hash_name => q{active_parameter}, } );
my $yaml_file_path   = catfile( $test_dir, q{ap_test.yaml} );

write_yaml(
    {
        yaml_href      => \%active_parameter,
        yaml_file_path => $yaml_file_path,
    }
);

## Then yaml file should exist
ok( -e $yaml_file_path, q{Created yaml file} );

open my $YAML, q{<}, $yaml_file_path
  or croak q{Cannot open}
  . $SPACE
  . $SINGLE_QUOTE
  . $DOUBLE_QUOTE
  . $yaml_file_path
  . $DOUBLE_QUOTE
  . $SINGLE_QUOTE
  . $COLON
  . $OS_ERROR
  . $NEWLINE;

## Given a number in hash
my $quoted_numeric_string;

LINE:
while ( my $line = <$YAML> ) {

    ($quoted_numeric_string) = $line =~ /(java_use_large_pages:\s+'1')/xsm;
    last LINE if ($quoted_numeric_string);
}

## Then number should be quoted in file
ok( $quoted_numeric_string, q{Found quoted numeric string} );

done_testing();
