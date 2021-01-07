package MIP::Get::File;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $COMMA $DOT $LOG_NAME $PIPE $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      get_io_files
      get_merged_infile_prefix
    };
}

sub get_io_files {

## Function : Get the io files per chain, id and stream
## Returns  : %io
## Arguments: $id             => Id (sample or case)
##          : $file_info_href => File info hash {REF}
##          : $parameter_href => Parameter hash {REF}
##          : $recipe_name    => Recipe name
##          : $stream         => Stream (in or out or temp)
##          : $temp_directory => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $id;
    my $file_info_href;
    my $parameter_href;
    my $recipe_name;
    my $stream;
    my $temp_directory;

    my $tmpl = {
        id => {
            defined     => 1,
            required    => 1,
            store       => \$id,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        stream => {
            allow       => [qw{ in temp out }],
            defined     => 1,
            required    => 1,
            store       => \$stream,
            strict_type => 1,
        },
        temp_directory => {
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Avoid autovivification of variable
    use Data::Diver qw{ Dive };
    use List::MoreUtils qw{ before };
    use MIP::Set::File qw{ set_io_files };

    ## Constants
    Readonly my $CHAIN_MAIN => q{CHAIN_MAIN};

    ## Unpack
    my $chain_id = $parameter_href->{$recipe_name}{chain};

    ## Not first in chain - return file features
    if ( Dive( $file_info_href, ( q{io}, $chain_id, $id, $recipe_name, $stream ) ) ) {

        return %{ $file_info_href->{io}{$chain_id}{$id}{$recipe_name} };
    }
    else {
        ## First in chain - need to find out stream file features of
        ## correct upstream recipe

        my $upstream_direction = q{out};

        ## Unpack
        my @order_recipes =
          @{ $parameter_href->{cache}{order_recipes_ref} };

        ## Find upstream recipes starting from (and not including) recipe_name
        my @upstream_recipes =
          reverse before { $_ eq $recipe_name } @order_recipes;

      UPSTREAM_RECIPE:
        foreach my $upstream_recipe (@upstream_recipes) {

            # Get chain id
            my $upstream_chain_id = $parameter_href->{$upstream_recipe}{chain};

            ## No io file features found in chain and stream
            next UPSTREAM_RECIPE
              if (
                not Dive(
                    $file_info_href,
                    ( q{io}, $upstream_chain_id, $id, $upstream_recipe, $upstream_direction )
                )
              );

            ## PARALLEL CHAIN with multiple recipes
            # second in chain
            if ( $upstream_chain_id eq $chain_id ) {

                ## Switch upstream out to recipe in - i.e. inherit from upstream
                _inherit_upstream_io_files(
                    {
                        chain_id           => $chain_id,
                        id                 => $id,
                        file_info_href     => $file_info_href,
                        recipe_name        => $recipe_name,
                        stream             => $stream,
                        temp_directory     => $temp_directory,
                        upstream_direction => $upstream_direction,
                        upstream_chain_id  => $upstream_chain_id,
                        upstream_recipe    => $upstream_recipe,
                    }
                );

                ##  Return set file features
                return %{ $file_info_href->{io}{$chain_id}{$id}{$recipe_name} };
            }

            ## Do not inherit from other chains than self or MAIN
            next UPSTREAM_RECIPE if ( $upstream_chain_id ne q{MAIN} );

            ## Found io file features found in chain, id, recipe and stream
            if (
                Dive(
                    $file_info_href,
                    ( q{io}, $upstream_chain_id, $id, $upstream_recipe, $upstream_direction )
                )
              )
            {

                ## Switch upstream out to recipe in - i.e. inherit from upstream
                _inherit_upstream_io_files(
                    {
                        chain_id           => $chain_id,
                        id                 => $id,
                        file_info_href     => $file_info_href,
                        recipe_name        => $recipe_name,
                        stream             => $stream,
                        temp_directory     => $temp_directory,
                        upstream_direction => $upstream_direction,
                        upstream_chain_id  => $upstream_chain_id,
                        upstream_recipe    => $upstream_recipe,
                    }
                );

                ##  Return set file features
                return %{ $file_info_href->{io}{$chain_id}{$id}{$recipe_name} };
            }
        }
    }

    ## At root of initation map - add base
    # Build infiles path
    my @base_file_paths =
      map { catfile( $file_info_href->{$id}{mip_infiles_dir}, $_ ) }
      @{ $file_info_href->{$id}{mip_infiles} };

    set_io_files(
        {
            chain_id       => $CHAIN_MAIN,
            id             => $id,
            file_paths_ref => \@base_file_paths,
            file_info_href => $file_info_href,
            recipe_name    => $recipe_name,
            stream         => $stream,
            temp_directory => $temp_directory,
        }
    );
    return %{ $file_info_href->{io}{$CHAIN_MAIN}{$id}{$recipe_name} };
}

sub get_merged_infile_prefix {

## Function : Get the merged infile prefix for sample id
## Returns  : $merged_infile_prefix
## Arguments: $file_info_href => File info hash {REF}
##          : $sample_id      => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $sample_id;

    my $tmpl = {
        file_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$file_info_href
        },
        sample_id => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$sample_id
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return $file_info_href->{$sample_id}{merged_infile};
}

sub _inherit_upstream_io_files {

## Function : Switch upstream out to recipe in - i.e. inherit from upstream
## Returns  : %io
## Arguments: $chain_id           => Chain id
##          : $id                 => Id (sample or case)
##          : $file_info_href     => File info hash {REF}
##          : $recipe_name        => Recipe name
##          : $stream             => Stream (in or out or temp)
##          : $temp_directory     => Temporary directory
##          : $upstream_direction => Upstream direction
##          : $upstream_chain_id  => Upstream chain id
##          : $upstream_recipe    => Upstream recipe

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $chain_id;
    my $id;
    my $file_info_href;
    my $recipe_name;
    my $stream;
    my $temp_directory;
    my $upstream_direction;
    my $upstream_chain_id;
    my $upstream_recipe;

    my $tmpl = {
        chain_id => {
            defined     => 1,
            required    => 1,
            store       => \$chain_id,
            strict_type => 1,
        },
        id => {
            defined     => 1,
            required    => 1,
            store       => \$id,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        stream => {
            allow       => [qw{ in temp out }],
            defined     => 1,
            required    => 1,
            store       => \$stream,
            strict_type => 1,
        },
        temp_directory => {
            store       => \$temp_directory,
            strict_type => 1,
        },
        upstream_direction => {
            defined     => 1,
            required    => 1,
            store       => \$upstream_direction,
            strict_type => 1,
        },
        upstream_chain_id => {
            defined     => 1,
            required    => 1,
            store       => \$upstream_chain_id,
            strict_type => 1,
        },
        upstream_recipe => {
            defined     => 1,
            required    => 1,
            store       => \$upstream_recipe,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Switch upstream out to recipe in - i.e. inherit from upstream
    my @upstream_outfile_paths =
      @{ $file_info_href->{io}{$upstream_chain_id}{$id}
          {$upstream_recipe}{$upstream_direction}{file_paths} };
    set_io_files(
        {
            chain_id       => $chain_id,
            id             => $id,
            file_paths_ref => \@upstream_outfile_paths,
            file_info_href => $file_info_href,
            recipe_name    => $recipe_name,
            stream         => $stream,
            temp_directory => $temp_directory,
        }
    );
    return;
}

1;
