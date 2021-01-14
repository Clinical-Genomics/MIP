package MIP::Validate::Data;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ %constraint };
}

our %constraint = (
    file_exists        => sub { return 1 if ( -e $_[0] );                  return; },
    file_is_executable => sub { return 1 if ( -x $_[0] );                  return; },
    dir_exists         => sub { return 1 if ( -d $_[0] );                  return; },
    is_digit           => sub { return 1 if ( $_[0] =~ / \A \d+ \z /sxm ); return; },
    is_gzipped         => sub { return 1 if ( $_[0] =~ / [.]gz$ /xms );    return 0; },
    plain_file_exists  => sub { return 1 if ( -f $_[0] );                  return; },
);

1;

