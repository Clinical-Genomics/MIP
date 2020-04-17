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

## MIPs lib/
use MIP::Constants
  qw{ $COLON $COMMA $DOUBLE_QUOTE $EQUALS $SEMICOLON @SINGULARITY_BIND_PATHS $WITH_SINGULARITY };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ parse_sing_bind_paths reduce_dir_paths };
}

sub parse_sing_bind_paths {

## Function : Parse singularity bind paths and add export command to array
## Returns  : $singularity_bind
## Arguments: $active_parameter_href       => The active parameters for this analysis hash {REF}
##          : $package_name                => Package name
##          : $source_environment_cmds_ref => Array with source environment commands {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $package_name;
    my $source_environment_cmds_ref;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        package_name => {
            defined     => 1,
            required    => 1,
            store       => \$package_name,
            strict_type => 1,
        },
        source_environment_cmds_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$source_environment_cmds_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return if not $WITH_SINGULARITY;

    my @export_bind_paths = @SINGULARITY_BIND_PATHS;

    ## Look for extra bind paths
    if ( $active_parameter_href->{singularity_recipe_bind_path}{$package_name} ) {

        ## Add extra paths
        push @export_bind_paths,
          @{ $active_parameter_href->{singularity_recipe_bind_path}{$package_name} };

        ## Check for redundant paths
        @export_bind_paths = reduce_dir_paths( { dir_paths_ref => \@export_bind_paths } );
    }

    my $singularity_bind = join $COMMA, @export_bind_paths;

    my $singularity_bind_var =
        q{export SINGULARITY_BIND}
      . $EQUALS
      . $DOUBLE_QUOTE
      . $singularity_bind
      . $DOUBLE_QUOTE
      . $SEMICOLON;

    push @{$source_environment_cmds_ref}, $singularity_bind_var;

    return;
}

sub reduce_dir_paths {

## Function : Parses directory paths and reduces them to a non-overlapping array. No check for existing files or directories
## Returns  : @reduced_dir_paths
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

    my @dir_paths;

    ## Split to dir path to array
  DIR_PATH:
    foreach my $dir_path ( @{$dir_paths_ref} ) {

        next DIR_PATH if ( not defined $dir_path );

        push @dir_paths, [ splitdir($dir_path) ];
    }

    ## Sort according to size
    @dir_paths = sort { @{$a} <=> @{$b} } @dir_paths;

    ## Reformat to strings
    @dir_paths = map { catdir( @{$_} ) } @dir_paths;

    my @reduced_dir_paths;

  BIND_PATH:
    while (@dir_paths) {

        ## Shift array
        my $dir_path = shift @dir_paths;

        ## Save path
        push @reduced_dir_paths, $dir_path;

        ## Get indexes of all the paths in the array that have an identical beginning to one we are testing
        ## The \Q and \E in the regex turns of interpolation
        my @match_idxs =
          grep { $dir_paths[$_] =~ / ^\Q$dir_path\E.* /xms } 0 .. $#dir_paths;

      MATCH_IDX:
        foreach my $match_idx ( reverse @match_idxs ) {

            ## Remove those paths with matching starts
            splice @dir_paths, $match_idx, 1;
        }
    }

    return @reduced_dir_paths;
}

1;
