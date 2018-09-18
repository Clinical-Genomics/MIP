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
use Modern::Perl qw{ 2014 };
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
Readonly my $FRW_SLASH => q{/};
Readonly my $COMMA     => q{,};
Readonly my $SPACE     => q{ };

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::Set::File}      => [qw{ set_io_files }],
        q{MIP::Test::Fixtures} => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::Set::File qw{ set_io_files };

diag(   q{Test set_io_files from File.pm v}
      . $MIP::Set::File::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Given
my $chain_id = 1;
my $id       = q{sample_1};
my @file_paths =
  ( catfile(qw{ a test dir file_1.txt}), catfile(qw{ a test dir file_2.txt}), );
my @file_path_prefixes =
  ( catfile(qw{ a test dir file_1}), catfile(qw{ a test dir file_2}), );
my @temp_file_paths =
  ( catfile(qw{ a temp dir file_1.txt}), catfile(qw{ a temp dir file_2.txt}), );
my @temp_file_path_prefixes =
  ( catfile(qw{ a temp dir file_1}), catfile(qw{ a temp dir file_2}), );
my %file_info;
my $stream         = q{in};
my $temp_directory = catfile(qw{ a temp dir});

set_io_files(
    {
        chain_id       => $chain_id,
        id             => $id,
        file_paths_ref => \@file_paths,
        file_info_href => \%file_info,
        stream         => $stream,
        temp_directory => $temp_directory,
    }
);

## Then set io features for $chain_id
is(
    $file_info{io}{$chain_id}{$id}{$stream}{dir_path},
    catfile(qw{ a test dir }) . $FRW_SLASH,
    q{Set file path}
);

is(
    $file_info{io}{$chain_id}{$id}{$stream}{dir_name},
    catfile(qw{ a test dir }),
    q{Set dir name}
);

is_deeply(
    $file_info{io}{$chain_id}{$id}{$stream}{file_names},
    [qw{ file_1.txt file_2.txt }],
    q{Set file name}
);

is_deeply( \@{ $file_info{io}{$chain_id}{$id}{$stream}{file_name_prefixes} },
    [qw{ file_1 file_2 }], q{Set file name prefixes} );

is_deeply( \@{ $file_info{io}{$chain_id}{$id}{$stream}{file_paths} },
    \@file_paths, q{Set file paths} );

is_deeply( \@{ $file_info{io}{$chain_id}{$id}{$stream}{file_path_prefixes} },
    \@file_path_prefixes, q{Set file path prefixes} );

is( $file_info{io}{$chain_id}{$id}{$stream}{file_suffix},
    q{.txt}, q{Set file suffix} );

## And for the temp dir
my $temp_stream = q{temp};

is(
    $file_info{io}{$chain_id}{$id}{$temp_stream}{dir_path},
    catfile(qw{ a temp dir }) . $FRW_SLASH,
    q{Set file path for temp stream}
);

is(
    $file_info{io}{$chain_id}{$id}{$temp_stream}{dir_name},
    catfile(qw{ a temp dir }),
    q{Set dir name for temp stream}
);

is_deeply(
    $file_info{io}{$chain_id}{$id}{$temp_stream}{file_names},
    [qw{ file_1.txt file_2.txt }],
    q{Set file name for temp stream}
);

is_deeply(
    \@{ $file_info{io}{$chain_id}{$id}{$temp_stream}{file_name_prefixes} },
    [qw{ file_1 file_2 }], q{Set file name prefixes for temp stream} );

is_deeply( \@{ $file_info{io}{$chain_id}{$id}{$temp_stream}{file_paths} },
    \@temp_file_paths, q{Set file paths for temp stream} );

is_deeply(
    \@{ $file_info{io}{$chain_id}{$id}{$temp_stream}{file_path_prefixes} },
    \@temp_file_path_prefixes, q{Set file path prefixes for temp stream} );

is( $file_info{io}{$chain_id}{$id}{$temp_stream}{file_suffix},
    q{.txt}, q{Set file suffix for temp stream} );
done_testing();
