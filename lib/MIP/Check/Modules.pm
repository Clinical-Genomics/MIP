package MIP::Check::Modules;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ parse_cpan_file };
}

sub parse_cpan_file {

## Function :
## Returns  : @cpanm_modules
## Arguments: $cpanfile_path => Path to cpanfile

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $cpanfile_path;

    my $tmpl = {
        cpanfile_path => {
            default     => 1,
            required    => 1,
            store       => \$cpanfile_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Module::CPANfile;
    use CPAN::Meta::Prereqs;

    ## Load cpanfile
    my $file = Module::CPANfile->load($cpanfile_path);
    ## Get hash_ref without objects
    my $file_href = $file->prereqs->as_string_hash;

    ## Get cpanm modules
    my @cpanm_modules = keys %{ $file_href->{runtime}{requires} };

    return @cpanm_modules;
}
1;
