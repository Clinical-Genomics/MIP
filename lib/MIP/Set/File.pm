package MIP::Set::File;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname fileparse };
use File::Spec::Functions qw{ catdir catfile splitpath };
use Params::Check qw{ check allow last_error };
use open qw{ :encoding(UTF-8) :std };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use List::MoreUtils qw { any uniq };
use Readonly;

BEGIN {

    use base qw{ Exporter };
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.08;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ set_file_compression_features set_file_prefix_tag set_infiles set_io_files set_merged_infile_prefix };
}

sub set_file_compression_features {

## Function : Set file compression features
## Returns  : $is_compressed, $read_file_command
## Arguments: $file_name => File name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_name;

    my $tmpl = {
        file_name => {
            defined     => 1,
            required    => 1,
            store       => \$file_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Path qw{ check_gzipped };

    my $read_file_command = q{gzip -d -c};

    ## Check if a file is gzipped.
    my $is_compressed = check_gzipped( { file_name => $file_name, } );

    ## Not compressed
    if ( not $is_compressed ) {

        ## File needs compression before starting analysis
        $read_file_command = q{cat};
    }
    return $is_compressed, $read_file_command;
}

sub set_file_prefix_tag {

## Function : Set the file tag depending on active recipes.
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $current_chain         => Name of current chain
##          : $case_id               => Family id {REF}
##          : $file_tag              => File tag to set
##          : $file_info_href        => Info on files hash {REF}
##          : $id                    => To change id for
##          : $recipe_name           => Recipe to add file tag for
##          : $temp_file_ending_href => Store sequential build of file tag

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $current_chain;
    my $file_tag;
    my $file_info_href;
    my $id;
    my $recipe_name;
    my $temp_file_ending_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        current_chain => {
            defined     => 1,
            required    => 1,
            store       => \$current_chain,
            strict_type => 1,
        },
        file_tag => {
            defined     => 1,
            required    => 1,
            store       => \$file_tag,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        id => {
            defined     => 1,
            required    => 1,
            store       => \$id,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        temp_file_ending_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$temp_file_ending_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## File_ending should be added for this recipe
    if ( $active_parameter_href->{$recipe_name} ) {

        _inherit_chain_main(
            {
                current_chain         => $current_chain,
                id                    => $id,
                temp_file_ending_href => $temp_file_ending_href,
            }
        );

        if ( defined $temp_file_ending_href->{$current_chain}{$id} ) {

            ## Add new file tag to build-up
            $file_info_href->{$id}{$recipe_name}{file_tag} =
              $temp_file_ending_href->{$current_chain}{$id} . $file_tag;
        }
        else {
            ## First recipe that should add filending

            $file_info_href->{$id}{$recipe_name}{file_tag} = $file_tag;
        }
    }
    else {
        ## Do not add recipe file_tag for this recipe but for previous recipes

        $file_info_href->{$id}{$recipe_name}{file_tag} =
          $temp_file_ending_href->{$current_chain}{$id};
    }

    return $file_info_href->{$id}{$recipe_name}{file_tag};
}

sub _inherit_chain_main {

## Function : Inherit file tags from MAIN chain.
## Returns  :
## Arguments: $current_chain         => Name of current chain
##          : $id                    => To change id for
##          : $temp_file_ending_href => Store sequential build of file tag

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $current_chain;
    my $id;
    my $temp_file_ending_href;

    my $tmpl = {
        current_chain => {
            defined     => 1,
            required    => 1,
            store       => \$current_chain,
            strict_type => 1,
        },
        id => {
            defined     => 1,
            required    => 1,
            store       => \$id,
            strict_type => 1,
        },
        temp_file_ending_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$temp_file_ending_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Check if branch and first recipe on branch
    if ( $current_chain ne q{MAIN}
        and not defined $temp_file_ending_href->{$current_chain}{$id} )
    {

        ## Inherit current MAIN chain.
        $temp_file_ending_href->{$current_chain}{$id} =
          $temp_file_ending_href->{MAIN}{$id};
    }
    return;
}

sub set_infiles {

## Function : Set the infile features i.e. dir and infiles
## Returns  :
## Arguments: $file_info_href   => File info hash {REF}
##          : $infile_directory => Infile directory
##          : $infiles_ref      => Infiles to check {REF}
##          : $sample_id        => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $infiles_ref;
    my $infile_directory;
    my $sample_id;

    my $tmpl = {
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        infiles_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$infiles_ref,
            strict_type => 1,
        },
        infile_directory => {
            defined     => 1,
            required    => 1,
            store       => \$infile_directory,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Set inputdir path hash
    $file_info_href->{$sample_id}{mip_infiles_dir} = $infile_directory;

    ## Set infiles in hash
    $file_info_href->{$sample_id}{mip_infiles} = [ @{$infiles_ref} ];
    return;
}

sub set_io_files {

## Function : Set the io files per chain and stream
## Returns  : io
## Arguments: $chain_id       => Chain of recipe
##          : $id             => Id (sample or case)
##          : $file_info_href => File info hash {REF}
##          : $file_paths_ref => File paths {REF}
##          : $recipe_name    => Recipe name
##          : $stream         => Stream (in or out or temp)
##          : $temp_directory => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $chain_id;
    my $id;
    my $file_info_href;
    my $file_paths_ref;
    my $recipe_name;
    my $stream;
    my $temp_directory;

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
        file_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$file_paths_ref,
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

    ## Alias
    my $io_recipe_href =
      \%{ $file_info_href->{io}{$chain_id}{$id}{$recipe_name} };

    ## Delete previous record (if any)
    delete $io_recipe_href->{$stream};

  FILE_PATH:
    foreach my $file_path ( @{$file_paths_ref} ) {

        my ( $file_name_prefix, $dirs, $suffix ) =
          fileparse( $file_path, qr/([.][^.]*)*/sxm );

        push @{ $io_recipe_href->{$stream}{file_names} },         basename($file_path);
        push @{ $io_recipe_href->{$stream}{file_name_prefixes} }, $file_name_prefix;
        push @{ $io_recipe_href->{$stream}{file_paths} },         $file_path;
        push @{ $io_recipe_href->{$stream}{file_path_prefixes} },
          catfile( $dirs, $file_name_prefix );

        ## Collect everything after first dot
        push @{ $io_recipe_href->{$stream}{file_suffixes} }, $suffix;
    }

    if ( scalar @{ $io_recipe_href->{$stream}{file_suffixes} } ) {

        _set_io_files_hash(
            {
                chain_id       => $chain_id,
                file_info_href => $file_info_href,
                id             => $id,
                recipe_name    => $recipe_name,
                stream         => $stream,
            }
        );

    }

    ## Split relative infile_path to file(s)
    my ( $infile_path_volume, $file_path_directory, $file_path_file_name ) =
      splitpath( $file_paths_ref->[0] );

    $io_recipe_href->{$stream}{dir_path} = $file_path_directory;
    $io_recipe_href->{$stream}{dir_path_prefix} =
      dirname( $file_paths_ref->[0] );

    ## Collect everything after last dot
    my ( $filename, $dirs, $suffix ) =
      fileparse( $file_paths_ref->[0], qr/[.][^.]*/sxm );
    $io_recipe_href->{$stream}{file_suffix} = $suffix;

    _set_io_files_constant(
        {
            chain_id       => $chain_id,
            file_info_href => $file_info_href,
            id             => $id,
            recipe_name    => $recipe_name,
            stream         => $stream,
        }
    );

    ## Also set the temporary file features for stream
    if ($temp_directory) {

        ## Switch to temp dir for path
        my @file_paths_temp =
          map { catfile( $temp_directory, $_ ) }
          @{ $io_recipe_href->{$stream}{file_names} };

        set_io_files(
            {
                chain_id       => $chain_id,
                id             => $id,
                file_paths_ref => \@file_paths_temp,
                file_info_href => $file_info_href,
                recipe_name    => $recipe_name,
                stream         => q{temp},
            }
        );
        return;
    }
    return;
}

sub set_merged_infile_prefix {

## Function : Set the merged infile prefix for sample id
## Returns  :
## Arguments: $file_info_href       => File info hash {REF}
##          : $merged_infile_prefix => Merged infile prefix
##          : $sample_id            => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $merged_infile_prefix;
    my $sample_id;

    my $tmpl = {
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        merged_infile_prefix => {
            defined     => 1,
            required    => 1,
            store       => \$merged_infile_prefix,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    $file_info_href->{$sample_id}{merged_infile} = $merged_infile_prefix;

    return;
}

sub _set_io_files_constant {

## Function : Set the io files per chain and stream for constant features
## Returns  : io
## Arguments: $chain_id       => Chain of recipe
##          : $id             => Id (sample or case)
##          : $file_info_href => File info hash {REF}
##          : $stream         => Stream (in or out or temp)

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $chain_id;
    my $id;
    my $file_info_href;
    my $recipe_name;
    my $stream;

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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Alias
    my $io_recipe_href =
      \%{ $file_info_href->{io}{$chain_id}{$id}{$recipe_name} };

    my %constant_map = (
        file_name_prefixes => q{file_name_prefix},
        file_paths         => q{file_path},
        file_path_prefixes => q{file_path_prefix},
        file_suffixes      => q{file_constant_suffix},
    );

    while ( my ( $file_feature, $file_constant_feature ) = each %constant_map ) {

        ## Get unique suffixes
        my @uniq_elements =
          uniq( @{ $io_recipe_href->{$stream}{$file_feature} } );

        ## If unique
        if ( scalar @uniq_elements == 1 ) {

            ## Set file constant suffix
            $io_recipe_href->{$stream}{$file_constant_feature} =
              $uniq_elements[0];
        }
    }
    return;
}

sub _set_io_files_hash {

## Function : Set the io hash files per chain, id, recipe and stream
## Returns  : io
## Arguments: $chain_id       => Chain of recipe
##          : $id             => Id (sample or case)
##          : $file_info_href => File info hash {REF}
##          : $file_paths_ref => File paths {REF}
##          : $recipe_name    => Recipe name
##          : $stream         => Stream (in or out or temp)

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $chain_id;
    my $id;
    my $file_info_href;
    my $file_paths_ref;
    my $recipe_name;
    my $stream;

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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Alias
    my $io_recipe_href =
      \%{ $file_info_href->{io}{$chain_id}{$id}{$recipe_name} };

    my %file_map = (
        file_names => q{file_name_href},
        file_paths => q{file_path_href},
    );

  FILE_MAP:
    while ( my ( $array_key, $hash_key ) = each %file_map ) {

      SUFFIX:
        while ( my ( $file_index, $suffix ) =
            each @{ $io_recipe_href->{$stream}{file_suffixes} } )
        {

            my $file = $io_recipe_href->{$stream}{$array_key}[$file_index];

            ## Find iterator string between dots
            my ($iterator) = $suffix =~ /[.]([^.]+)[.]/sxm;

            if ($iterator) {

                $io_recipe_href->{$stream}{$hash_key}{$iterator} = $file;
            }
        }
    }
    return;
}

1;
