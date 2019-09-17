package MIP::Parse::Singularity;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir splitdir };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ parse_sing_bind_paths };
}

sub parse_sing_bind_paths {

## Function : Parse singularity bind paths and reduces them to a non-overlapping array
## Returns  : @reduced_bind_paths
## Arguments: $dir_paths_ref => Directory paths to parse {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $dir_paths_ref;

    my $tmpl = {
        dir_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$dir_paths_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @bind_paths;

    ## Split to dir path to array
    foreach my $dir_path ( @{$dir_paths_ref} ) {

        push @bind_paths, [ splitdir($dir_path) ];
    }

    ## Sort according to size
    @bind_paths = sort { @{$a} <=> @{$b} } @bind_paths;

    ## Reformat to strings
    @bind_paths = map { catdir( @{$_} ) } @bind_paths;

    my @reduced_bind_paths;

  BIND_PATH:
    while (@bind_paths) {

        ## Shift array
        my $bind_path = shift @bind_paths;

        ## Save path
        push @reduced_bind_paths, $bind_path;

        ## Get indexes of all the paths in the array that have an identical beginning to one we are testing
        ## The \Q and \E in the regex turns of interpolation
        my @match_idxs =
          grep { $bind_paths[$_] =~ / ^\Q$bind_path\E.* /xms } 0 .. $#bind_paths;

      MATCH_IDX:
        foreach my $match_idx ( reverse @match_idxs ) {

            ## Remove those paths with matching starts
            splice @bind_paths, $match_idx, 1;
        }
    }

    return @reduced_bind_paths;
}

1;
