package MIP::Parse::Singularity;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile splitdir };
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
use MIP::Constants qw{ $SPACE @SINGULARITY_BIND_PATHS %SINGULARITY_CONTAINER };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ parse_commands_for_singularity parse_sing_bind_paths };
}

sub parse_commands_for_singularity {

## Function : Parse command array from unix_write_to_file
## Returns  :
## Arguments: $commands_ref => Commands to write to file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $commands_ref;

    my $tmpl = {
        commands_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$commands_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Program::Singularity qw{ singularity_exec };

    # Return if empty command array
    return if ( not @{$commands_ref} );

    # Get possible executable to test
    my $program_executable = $commands_ref->[0];

    ## Check if executable has a singularity container
    my $singularity_container;

    if ( defined $SINGULARITY_CONTAINER{$program_executable} ) {

        $singularity_container =
          $SINGULARITY_CONTAINER{$program_executable}{container_path};
    }

    # Return if no singularity image was found
    return if ( not $singularity_container );

    ## Get bind paths
    my @bind_paths = @SINGULARITY_BIND_PATHS;
    if ( $SINGULARITY_CONTAINER{$program_executable}{extra_bind_paths} ) {

        push @bind_paths,
          @{ $SINGULARITY_CONTAINER{$program_executable}{extra_bind_paths} };
    }

    ## Remove potential undef values
    @bind_paths = grep { defined } @bind_paths;

    ## Remove overlapping paths
    @bind_paths = parse_sing_bind_paths( { dir_paths_ref => \@bind_paths, } );

    # Build complete command
    @{$commands_ref} = singularity_exec(
        {
            bind_paths_ref                 => \@bind_paths,
            singularity_container          => $singularity_container,
            singularity_container_cmds_ref => $commands_ref,
        }
    );

    return @{$commands_ref};
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
