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
use MIP::Constants qw{ $COLON $COMMA $SPACE $TAB };
use MIP::Test::Fixtures qw{ test_standard_cli };

my $VERBOSE = 1;
our $VERSION = 1.02;

$VERBOSE = test_standard_cli(
    {
        verbose => $VERBOSE,
        version => $VERSION,
    }
);

## Constants
Readonly my $INDEX_READ_1_SECOND_PAIR => 2;
Readonly my $INDEX_READ_2_SECOND_PAIR => 3;

BEGIN {

    use MIP::Test::Fixtures qw{ test_import };

### Check all internal dependency modules and imports
## Modules with import
    my %perl_module = (
        q{MIP::File::Format::Star_fusion} => [qw{ create_star_fusion_sample_file }],
        q{MIP::Test::Fixtures}            => [qw{ test_standard_cli }],
    );

    test_import( { perl_module_href => \%perl_module, } );
}

use MIP::File::Format::Star_fusion qw{ create_star_fusion_sample_file };

diag(   q{Test create_star_fusion_sample_file from Star_fusion.pm v}
      . $MIP::File::Format::Star_fusion::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

my $test_dir = File::Temp->newdir();

# Create anonymous filehandle
my $filehandle = IO::Handle->new();

# For storing info to write
my $file_content;

## Given a sample id and infiles
my @infile_paths = (
    catfile(qw{ a dir file_1.fastq }),   catfile(qw{ a dir file_2.fastq }),
    catfile(qw{ a dir file_x_1.fastq }), catfile(qw{ a dir file_x_2.fastq }),
);
my $sample_id         = q{sample_1};
my $samples_file_path = catfile( $test_dir, q{sample_file} );
my %file_info         = (
    $sample_id => {
        file_prefix_no_direction => {
            file   => q{paired-end},
            file_x => q{paired-end}
        },
    },
);

## Store file content in memory by using referenced variable
open $filehandle, q{>}, \$file_content
  or croak q{Cannot write to} . $SPACE . $file_content . $COLON . $SPACE . $OS_ERROR;

create_star_fusion_sample_file(
    {
        filehandle        => $filehandle,
        file_info_href    => \%file_info,
        infile_paths_ref  => \@infile_paths,
        samples_file_path => $samples_file_path,
        sample_id         => $sample_id,
    }
);
close $filehandle;

## Then file content should exist and string with sample_id and path should written
ok( $file_content, q{Created file content} );

my ($returned_sample_id) = $file_content =~ /($sample_id)/msx;

my ($returned_samples_file_path) = $file_content =~ /($samples_file_path)/mxs;

is( $returned_sample_id, $sample_id, q{Found sample_id in line} );

is( $returned_samples_file_path, $samples_file_path, q{Found sample path in line} );
my $second_pair = $infile_paths[$INDEX_READ_1_SECOND_PAIR] . q{\S+}
  . $infile_paths[$INDEX_READ_2_SECOND_PAIR];
my ($returned_second_pair) = $file_content =~ /($second_pair)/msx;

ok( $returned_second_pair, q{Found second pair} );

done_testing();
