package MIP::File_info;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname fileparse };
use File::Spec::Functions qw{ catdir catfile splitpath };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use List::MoreUtils qw { uniq };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $DOT $EMPTY_STR $GENOME_VERSION $LOG_NAME $NEWLINE $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      add_sample_fastq_file_lanes
      add_sample_no_direction_infile_prefixes
      check_parameter_metafiles
      get_consensus_sequence_run_type
      get_io_files
      get_is_sample_files_compressed
      get_merged_infile_prefix
      get_sample_file_attribute
      get_sample_fastq_file_lanes
      get_sampling_fastq_files
      parse_file_compression_features
      parse_files_compression_status
      parse_io_outfiles
      parse_sample_fastq_file_attributes
      parse_select_file_contigs
      set_alt_loci_contigs
      set_dict_contigs
      set_file_tag
      set_human_genome_reference_features
      set_infiles
      set_io_files
      set_is_sample_files_compressed
      set_human_genome_reference_features
      set_merged_infile_prefix
      set_primary_contigs
      set_sample_file_attribute
      set_sample_max_parallel_processes_count
      set_select_file_contigs
    };
}

## Constants
Readonly my $INTERLEAVED_READ_DIRECTION => 3;

sub add_sample_fastq_file_lanes {

## Function : Add sample fastq file lane to lanes
## Returns  :
## Arguments: $direction      => Read direction
##          : $file_info_href => File info hash {REF}
##          : $lane           => Lane number
##          : $sample_id      => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $direction;
    my $file_info_href;
    my $lane;
    my $sample_id;

    my $tmpl = {
        direction => {
            allow       => [ undef, 1, 2, $INTERLEAVED_READ_DIRECTION, ],
            required    => 1,
            store       => \$direction,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        lane => {
            required    => 1,
            store       => \$lane,
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

    return if ( not $lane );

    return if ( not $direction == 1 );

    ## Add lane
    push @{ $file_info_href->{$sample_id}{lanes} }, $lane;

    return;
}

sub add_sample_no_direction_infile_prefixes {

## Function : Add sample fastq file prefix without read direction in file name
## Returns  :
## Arguments: $file_info_href  => File info hash {REF}
##          : $mip_file_format => Mip file format without read direction and ".fastq(.gz)"
##          : $sample_id       => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $mip_file_format;
    my $sample_id;

    my $tmpl = {
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        mip_file_format => {
            defined     => 1,
            required    => 1,
            store       => \$mip_file_format,
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

    ## Add no_direction_infile_prefixes
    push @{ $file_info_href->{$sample_id}{no_direction_infile_prefixes} }, $mip_file_format;

    return;
}

sub check_parameter_metafiles {

## Function : Checks parameter metafile exists
## Returns  :
## Arguments: $active_parameter_href => Holds all set parameter for analysis
##          : $file_info_href        => File info hash {REF}
##          : $parameter_href        => Holds all parameters

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Reference qw{ parse_meta_file_suffixes };
    use MIP::Parameter qw{ get_parameter_attribute };

  PARAMETER:
    foreach my $parameter_name ( keys %{$file_info_href} ) {

        ## Active parameter
        my $parameter = $active_parameter_href->{$parameter_name};

        next PARAMETER if ( not $parameter );

        my @associated_recipes = get_parameter_attribute(
            {
                attribute      => q{associated_recipe},
                parameter_href => $parameter_href,
                parameter_name => $parameter_name,
            }
        );
        ## Find any active recipe among associated recipes
        my $has_active_recipe =
          grep { defined and $_ >= 1 } @{$active_parameter_href}{@associated_recipes};

        next PARAMETER if ( not $has_active_recipe );

        if ( ref $parameter eq q{HASH} ) {

          PATH:
            for my $path ( keys %{$parameter} ) {

                ## Checks files to be built by combining filename stub with fileendings
                parse_meta_file_suffixes(
                    {
                        active_parameter_href  => $active_parameter_href,
                        file_name              => $path,
                        meta_file_suffixes_ref => \@{ $file_info_href->{$parameter_name} },
                        parameter_href         => $parameter_href,
                        parameter_name         => $parameter_name,
                    }
                );

                ## If single $path needs building - build for all as switch
                ## is set on parameter_name and not path
                my $build_status = get_parameter_attribute(
                    {
                        attribute      => q{build_file},
                        parameter_href => $parameter_href,
                        parameter_name => $parameter_name,
                    }
                );
                next PARAMETER if ($build_status);
            }
            next PARAMETER;
        }

        ## Checks files to be built by combining filename stub with fileendings
        parse_meta_file_suffixes(
            {
                active_parameter_href  => $active_parameter_href,
                file_name              => $parameter,
                meta_file_suffixes_ref => \@{ $file_info_href->{$parameter_name} },
                parameter_href         => $parameter_href,
                parameter_name         => $parameter_name,
            }
        );
    }
    return;
}

sub get_consensus_sequence_run_type {

## Function : Get consensus sequence run type across samples
## Returns  : 0 | $consensus_type
## Arguments: $file_info_href  => File info hash {REF}
##          : $sample_ids_ref  => Sample ids

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $sample_ids_ref;

    my $tmpl = {
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        sample_ids_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$sample_ids_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $has_consensus;
    my $consensus_type;

    ## Get sequence run modes
  SAMPLE_ID:
    foreach my $sample_id ( @{$sample_ids_ref} ) {

        my %seen;

        my %file_info_sample = get_sample_file_attribute(
            {
                file_info_href => $file_info_href,
                sample_id      => $sample_id,
            }
        );

      INFILE_PREFIX:
        foreach my $infile_prefix ( @{ $file_info_sample{no_direction_infile_prefixes} } ) {

            my $sequence_run_type = get_sample_file_attribute(
                {
                    attribute      => q{sequence_run_type},
                    file_info_href => $file_info_href,
                    file_name      => $infile_prefix,
                    sample_id      => $sample_id,
                }
            );
            $seen{$sequence_run_type} = $sequence_run_type;
            $consensus_type = $sequence_run_type;
        }

        ## Turn of recipe if multiple sequence run types are present
        $has_consensus = scalar keys %seen <= 1 ? 1 : 0;
        return 0 if ( not $has_consensus );
    }
    return $consensus_type;
}

sub get_io_files {

## Function : Get the io files per chain, id and stream
## Returns  : %io
## Arguments: $file_info_href => File info hash {REF}
##          : $id             => Id (sample or case)
##          : $parameter_href => Parameter hash {REF}
##          : $recipe_name    => Recipe name
##          : $stream         => Stream (in or out or temp)
##          : $temp_directory => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $id;
    my $parameter_href;
    my $recipe_name;
    my $stream;
    my $temp_directory;

    my $tmpl = {
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
            allow       => [qw{ in out temp }],
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

    use Data::Diver qw{ Dive };
    use List::MoreUtils qw{ before };
    use MIP::File_info qw{ set_io_files };
    use MIP::Parameter qw{ get_parameter_attribute };

    ## Constants
    Readonly my $CHAIN_MAIN         => q{CHAIN_MAIN};
    Readonly my $UPSTREAM_DIRECTION => q{out};

    ## Unpack
    my $chain_id = get_parameter_attribute(
        {
            attribute      => q{chain},
            parameter_href => $parameter_href,
            parameter_name => $recipe_name,
        }
    );

    ## Not first in chain - return file features
    if ( Dive( $file_info_href, ( q{io}, $chain_id, $id, $recipe_name, $stream ) ) ) {

        return %{ $file_info_href->{io}{$chain_id}{$id}{$recipe_name} };
    }
    else {
        ## First in chain - need to find out stream file features of correct upstream recipe

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
                    ( q{io}, $upstream_chain_id, $id, $upstream_recipe, $UPSTREAM_DIRECTION )
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
                        upstream_direction => $UPSTREAM_DIRECTION,
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
                    ( q{io}, $upstream_chain_id, $id, $upstream_recipe, $UPSTREAM_DIRECTION )
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
                        upstream_direction => $UPSTREAM_DIRECTION,
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

sub get_is_sample_files_compressed {

## Function : Get sample files compression status
## Returns  : 0 | 1
## Arguments: $file_info_href  => File info hash {REF}
##          : $sample_id       => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $sample_id;

    my $tmpl = {
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
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

    use Data::Diver qw{ Dive };

    if ( defined Dive( $file_info_href, ( q{is_files_compressed}, $sample_id ) ) ) {

        ## Return files compression status
        return $file_info_href->{is_files_compressed}{$sample_id};
    }
    return;
}

sub get_merged_infile_prefix {

## Function : Get the merged infile prefix for a sample id
## Returns  : $merged_infile_prefix
## Arguments: $file_info_href => File info hash {REF}
##          : $sample_id      => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $sample_id;

    my $tmpl = {
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
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

    return $file_info_href->{$sample_id}{merged_infile};
}

sub get_sample_file_attribute {

## Function : Get sample file attributes
## Returns  : %{ $file_info_href->{$sample_id} } | %{ $file_info_href->{$sample_id}{$file_name} } | $attribute
## Arguments: $attribute       => Attribute key
##          : $file_info_href  => File info hash {REF}
##          : $file_name       => File name
##          : $sample_id       => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $attribute;
    my $file_info_href;
    my $file_name;
    my $sample_id;

    my $tmpl = {
        attribute => {
            allow => [
                qw{ date
                  direction
                  flowcell
                  index
                  infile_sample_id
                  is_file_compressed
                  is_interleaved
                  lane
                  read_file_command
                  read_length
                  sequence_run_type
                  }
            ],
            store       => \$attribute,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        file_name => {
            store       => \$file_name,
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

    if ( not $attribute and not $file_name ) {

        ## Return entire sample_id hash
        return %{ $file_info_href->{$sample_id} };
    }
    if ( not $attribute ) {

        ## Return entire file name array
        if ( ref $file_info_href->{$sample_id}{$file_name} eq q{ARRAY} ) {
            return @{ $file_info_href->{$sample_id}{$file_name} };
        }

        ## Return entire file name hash
        if ( ref $file_info_href->{$sample_id}{$file_name} eq q{HASH} ) {
            return %{ $file_info_href->{$sample_id}{$file_name} };
        }

    }
    ## Get attribute
    my $stored_attribute =
      $file_info_href->{$sample_id}{$file_name}{$attribute};

    ## Return requested attribute
    return $stored_attribute;
}

sub get_sample_fastq_file_lanes {

## Function : Get sample fastq filelanes
## Returns  :
## Arguments: $file_info_href => File info hash {REF}
##          : $sample_id      => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $sample_id;

    my $tmpl = {
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
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

    return @{ $file_info_href->{$sample_id}{lanes} };
}

sub get_sampling_fastq_files {

## Function : Get sample fastq files. Either single-end, paired-end or interleaved
## Returns  : $is_interleaved_fastq, @fastq_files
## Arguments: $file_info_sample_href => File info sample hash

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_sample_href;

    my $tmpl = {
        file_info_sample_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_sample_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @fastq_files;

    ## Perform per single-end or read pair
    my $paired_end_tracker = 0;

  INFILE_PREFIX:
    foreach my $infile_prefix ( @{ $file_info_sample_href->{no_direction_infile_prefixes} } ) {

        push @fastq_files, $file_info_sample_href->{mip_infiles}[$paired_end_tracker];

        my $sequence_run_type = $file_info_sample_href->{$infile_prefix}{sequence_run_type};

        # If second read direction is present
        if ( $sequence_run_type eq q{paired-end} ) {

            # Increment to collect correct read 2
            $paired_end_tracker = $paired_end_tracker + 1;
            push @fastq_files, $file_info_sample_href->{mip_infiles}[$paired_end_tracker];
        }

        my $is_interleaved_fastq = $sequence_run_type eq q{interleaved} ? 1 : 0;

        ## Only perform once per sample and fastq file(s)
        return $is_interleaved_fastq, @fastq_files;
    }
    return;
}

sub parse_file_compression_features {

## Function : Parse file compression features
## Returns  : $attribute{read_file_command}
## Arguments: $file_info_href => File info hash {REF}
##          : $file_name      => File name
##          : $sample_id      => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $file_name;
    my $sample_id;

    my $tmpl = {
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        file_name => {
            defined     => 1,
            required    => 1,
            store       => \$file_name,
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

    use MIP::Validate::Data qw{ %constraint };

    my %attribute = (
        is_file_compressed => 0,
        read_file_command  => q{cat},
    );

    if ( $constraint{is_gzipped}->($file_name) ) {

        $attribute{is_file_compressed} = 1;
        $attribute{read_file_command}  = q{gzip -d -c};
    }

  ATTRIBUTES:
    while ( my ( $attribute, $attribute_value ) = each %attribute ) {

        set_sample_file_attribute(
            {
                attribute       => $attribute,
                attribute_value => $attribute_value,
                file_info_href  => $file_info_href,
                file_name       => $file_name,
                sample_id       => $sample_id,
            }
        );
    }
    return $attribute{read_file_command};
}

sub parse_files_compression_status {

## Function : Parse files compression status
## Returns  :
## Arguments: $file_info_href => File info hash {REF}
##          : $sample_id      => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $sample_id;

    my $tmpl = {
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
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

    my $is_compressed      = 0;
    my $compression_status = 0;

    ## Unpack
    my @infiles = @{ $file_info_href->{$sample_id}{mip_infiles} };

  FILE_NAME:
    foreach my $file_name (@infiles) {

        my $is_file_compressed = get_sample_file_attribute(
            {
                attribute      => q{is_file_compressed},
                file_info_href => $file_info_href,
                file_name      => $file_name,
                sample_id      => $sample_id,
            }
        );
        next FILE_NAME if ( not $is_file_compressed );

        $is_compressed++;
    }
    if ( $is_compressed == @infiles ) {
        $compression_status = 1;
    }
    ## Set is_files_compressed per sample global boolean
    set_is_sample_files_compressed(
        {
            compression_status => $compression_status,
            file_info_href     => $file_info_href,
            sample_id          => $sample_id,
        }
    );
    return;
}

sub parse_io_outfiles {

## Function : Set and get the io files per chain, id and stream
## Returns  : %io
## Arguments: $chain_id               => Chain of recipe
##          : $file_info_href         => File info hash {REF}
##          : $file_name_prefixes     => Build outfile using file name prefix
##          : $file_name_prefixes_ref => Build outfile using file name prefixes {REF}
##          : $file_paths_ref         => File paths {REF}
##          : $id                     => Id (sample or case)
##          : $iterators_ref          => Build outfile using iterator (e.g contigs) {REF}
##          : $outdata_dir            => Outdata directory
##          : $parameter_href         => Parameter hash {REF}
##          : $recipe_name            => Recipe name
##          : $stream                 => Stream (out)
##          : $temp_directory         => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $chain_id;
    my $file_info_href;
    my $file_name_prefix;
    my $file_name_prefixes_ref;
    my $file_paths_ref;
    my $id;
    my $iterators_ref;
    my $outdata_dir;
    my $parameter_href;
    my $recipe_name;
    my $temp_directory;

    ## Default(s)
    my $stream;

    my $tmpl = {
        chain_id => {
            defined     => 1,
            required    => 1,
            store       => \$chain_id,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        file_name_prefix => {
            store       => \$file_name_prefix,
            strict_type => 1,
        },
        file_name_prefixes_ref => {
            default     => [],
            store       => \$file_name_prefixes_ref,
            strict_type => 1,
        },
        file_paths_ref => {
            default     => [],
            store       => \$file_paths_ref,
            strict_type => 1,
        },
        id => {
            defined     => 1,
            required    => 1,
            store       => \$id,
            strict_type => 1,
        },
        iterators_ref => {
            default     => [],
            store       => \$iterators_ref,
            strict_type => 1,
        },
        outdata_dir => {
            store       => \$outdata_dir,
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
            allow       => [qw{ out }],
            default     => q{out},
            store       => \$stream,
            strict_type => 1,
        },
        temp_directory => {
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File_info qw{ get_io_files set_io_files };
    use MIP::Parameter qw{ get_recipe_attributes };

    ## Build default @file_paths
    if ( not @{$file_paths_ref} and $outdata_dir ) {

        my %recipe = get_recipe_attributes(
            {
                parameter_href => $parameter_href,
                recipe_name    => $recipe_name,
            }
        );
        my $outfile_tag    = $recipe{file_tag}       //= $EMPTY_STR;
        my $outfile_suffix = $recipe{outfile_suffix} //= $EMPTY_STR;
        my $directory      = catdir( $outdata_dir, $id, $recipe_name );

        ## Default paths with iterators
        if ( @{$iterators_ref} and $file_name_prefix ) {

            ## Localize as we will mutate elements
            my @iterators = @{$iterators_ref};
            foreach my $iterator (@iterators) {

                ## Add "." if not empty string
                $iterator = $iterator ne $EMPTY_STR ? $DOT . $iterator : $iterator;
            }
            @{$file_paths_ref} =
              map { catfile( $directory, $file_name_prefix . $outfile_tag . $_ . $outfile_suffix ) }
              @iterators;
        }
        ## Default paths without iterators
        else {

            ## $file_name_prefixes_needs to be set
            croak q{Missing argument!} if not @{$file_name_prefixes_ref};
            @{$file_paths_ref} =
              map { catfile( $directory, $_ . $outfile_tag . $outfile_suffix ) }
              @{$file_name_prefixes_ref};
        }
    }

    ## Set the io files per chain and stream
    set_io_files(
        {
            chain_id       => $chain_id,
            file_paths_ref => $file_paths_ref,
            file_info_href => $file_info_href,
            id             => $id,
            recipe_name    => $recipe_name,
            stream         => $stream,
            temp_directory => $temp_directory,
        }
    );

    my %io = get_io_files(
        {
            file_info_href => $file_info_href,
            id             => $id,
            parameter_href => $parameter_href,
            recipe_name    => $recipe_name,
            stream         => $stream,
        }
    );

    return %io;
}

sub parse_sample_fastq_file_attributes {

## Function : Parse sample fastq file attributes
## Returns  : %infile_info
## Arguments: $file_info_href => File info hash {REF}
##          : $file_name      => Fast file name
##          : $infiles_dir    => Sample infile dir of fastq files
##          : $sample_id      => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $file_name;
    my $infiles_dir;
    my $sample_id;

    my $tmpl = {
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        file_name => {
            defined     => 1,
            required    => 1,
            store       => \$file_name,
            strict_type => 1,
        },
        infiles_dir => {
            defined     => 1,
            required    => 1,
            store       => \$infiles_dir,
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

    use MIP::Fastq qw{ check_interleaved get_read_length parse_fastq_infiles_format };

    ## Parse infile according to filename convention
    my %infile_info = parse_fastq_infiles_format(
        {
            file_name => $file_name,
            sample_id => $sample_id,
        }
    );

    ## Parse compression features
    $infile_info{read_file_command} = parse_file_compression_features(
        {
            file_info_href => $file_info_href,
            file_name      => $file_name,
            sample_id      => $sample_id,
        }
    );

    ## Get sequence read length from file
    $infile_info{read_length} = get_read_length(
        {
            file_path         => catfile( $infiles_dir, $file_name ),
            read_file_command => $infile_info{read_file_command},
        }
    );

    ## Is file interleaved and have proper read direction
    $infile_info{is_interleaved} = check_interleaved(
        {
            file_path         => catfile( $infiles_dir, $file_name ),
            read_file_command => $infile_info{read_file_command},
        }
    );

    add_sample_fastq_file_lanes(
        {
            direction      => $infile_info{direction},
            file_info_href => $file_info_href,
            lane           => $infile_info{lane},
            sample_id      => $sample_id,
        }
    );

    ## Transfer to file_info hash
  ATTRIBUTE:
    while ( my ( $attribute, $attribute_value ) = each %infile_info ) {

        set_sample_file_attribute(
            {
                attribute       => $attribute,
                attribute_value => $attribute_value,
                file_info_href  => $file_info_href,
                file_name       => $file_name,
                sample_id       => $sample_id,
            }
        );
    }
    return %infile_info;
}

sub parse_select_file_contigs {

## Function : Parse select file contigs
## Returns  :
## Arguments: $consensus_analysis_type => Consensus analysis type for checking e.g. WGS specific files
##          : $file_info_href          => File info hash {REF}
##          : $select_file_path        => Select file path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $consensus_analysis_type;
    my $file_info_href;
    my $select_file_path;

    my $tmpl = {
        consensus_analysis_type => {
            defined     => 1,
            required    => 1,
            store       => \$consensus_analysis_type,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        select_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$select_file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Contigs qw{ check_select_file_contigs sort_contigs_to_contig_set };
    use MIP::Reference qw{ get_select_file_contigs };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    if ($select_file_path) {

        ## Collects sequences contigs used in select file
        my @select_file_contigs = get_select_file_contigs(
            {
                select_file_path => $select_file_path,
            }
        );

        ## Set in file_info hash
        set_select_file_contigs(
            {
                file_info_href          => $file_info_href,
                select_file_contigs_ref => \@select_file_contigs,
            }
        );
        ## Check that select file contigs is a subset of primary contigs
        check_select_file_contigs(
            {
                contigs_ref             => $file_info_href->{contigs},
                select_file_contigs_ref => $file_info_href->{select_file_contigs},
            }
        );

        ## Sorts array depending on reference array. NOTE: Only entries present in reference array will survive in sorted array.
        my %contig_sort_map = (
            select_file_contigs        => q{contigs},
            sorted_select_file_contigs => q{contigs_size_ordered},
        );
        while ( my ( $contigs_set_name, $sort_reference ) = each %contig_sort_map ) {

            @{ $file_info_href->{$contigs_set_name} } = sort_contigs_to_contig_set(
                {
                    consensus_analysis_type    => $consensus_analysis_type,
                    sort_contigs_ref           => $file_info_href->{select_file_contigs},
                    sort_reference_contigs_ref => $file_info_href->{$sort_reference},
                }
            );
        }
    }
    return 1;
}

sub set_alt_loci_contigs {

## Function : Set alternative loci contigs
## Returns  :
## Arguments: $alt_contig_set_name => Alt contig set identifier
##          : $file_info_href      => File info hash {REF}
##          : $primary_contig_href => Primary contig hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $alt_contig_set_name;
    my $file_info_href;
    my $primary_contig_href;

    my $tmpl = {
        alt_contig_set_name => {
            defined     => 1,
            required    => 1,
            store       => \$alt_contig_set_name,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        primary_contig_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$primary_contig_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Set alternative loci contig set
    @{ $file_info_href->{$alt_contig_set_name} } =
      grep { not exists $primary_contig_href->{$_} } @{ $file_info_href->{dict_contigs} };
    return;
}

sub set_primary_contigs {

## Function : Set primary contigs
## Returns  :
## Arguments: $contig_set_name     => Contig set identifier
##          : $file_info_href      => File info hash {REF}
##          : $primary_contigs_ref => Primary contig hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $contig_set_name;
    my $file_info_href;
    my $primary_contigs_ref;

    my $tmpl = {
        contig_set_name => {
            defined     => 1,
            required    => 1,
            store       => \$contig_set_name,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        primary_contigs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$primary_contigs_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Set bam contig sets
    @{ $file_info_href->{$contig_set_name} } = @{$primary_contigs_ref};

    return;
}

sub set_dict_contigs {

## Function : Set sequence contigs used in analysis from human genome sequence
##          : dictionnary (.dict file)
## Returns  :
## Arguments: $dict_file_path => Dict file path
##          : $file_info_href => File info hash {REF}
##          : $parameter_href => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $dict_file_path;
    my $file_info_href;
    my $parameter_href;

    my $tmpl = {
        dict_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$dict_file_path,
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Parameter qw{ get_parameter_attribute };
    use MIP::Reference qw{ get_dict_contigs };

    ## Get sequence contigs from human reference ".dict" file since it exists
    my $build_status = get_parameter_attribute(
        {
            attribute      => q{build_file},
            parameter_href => $parameter_href,
            parameter_name => q{human_genome_reference_file_endings},
        }
    );

## File needs to be built before getting contigs
    return if ($build_status);

    @{ $file_info_href->{dict_contigs} } = get_dict_contigs(
        {
            dict_file_path => $dict_file_path,
        }
    );

    return;
}

sub set_file_tag {

## Function : Set the file tag depending on id, branch and recipe
## Returns  :
## Arguments: $file_info_href => Info on files hash {REF}
##          : $file_tag       => File tag to set
##          : $id             => To change id for case or sample
##          : $recipe_name    => Recipe to add file tag for

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $file_tag;
    my $id;
    my $recipe_name;

    my $tmpl = {
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        file_tag => {
            required    => 1,
            store       => \$file_tag,
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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    $file_info_href->{$id}{$recipe_name}{file_tag} = $file_tag;

    return;
}

sub set_human_genome_reference_features {

## Function : Detect version and source of the human_genome_reference: Source (hg19 or grch) as well as compression status.
##            Used to change capture kit genome reference version later
## Returns  :
##          : $file_info_href         => File info hash {REF}
##          : $human_genome_reference => The human genome
##          : $parameter_href         => Parameter hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $human_genome_reference;
    my $parameter_href;

    my $tmpl = {
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        human_genome_reference => {
            defined     => 1,
            required    => 1,
            store       => \$human_genome_reference,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Constants qw{ set_genome_build_constants };
    use MIP::Parameter qw{ set_parameter_build_file_status };
    use MIP::Validate::Data qw{ %constraint };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Different regexes for the two sources.
    ## i.e. Don't allow subversion of Refseq genome
    my %genome_source = (
        grch => qr/grch(\d+[.]\d+ | \d+)/xsm,
        hg   => qr/hg(\d+)/xsm,
    );

  GENOME_PREFIX:
    foreach my $genome_prefix ( keys %genome_source ) {

        ## Capture version
        my ($genome_version) =
          $human_genome_reference =~ m/ $genome_source{$genome_prefix}_homo_sapiens /xms;

        next GENOME_PREFIX if ( not $genome_version );

        $file_info_href->{human_genome_reference_version} = $genome_version;
        $file_info_href->{human_genome_reference_source}  = $genome_prefix;

        ## Only set global constant once
        last GENOME_PREFIX if ($GENOME_VERSION);

        set_genome_build_constants(
            {
                genome_version => $genome_version,
                genome_source  => $genome_prefix,
            }
        );
        last GENOME_PREFIX;
    }
    if ( not $file_info_href->{human_genome_reference_version} ) {

        $log->fatal(
                q{MIP cannot detect what version of human_genome_reference you have supplied.}
              . $SPACE
              . q{Please supply the reference on this format: [sourceversion]_[species] e.g. 'grch37_homo_sapiens' or 'hg19_homo_sapiens'}
              . $NEWLINE );
        exit 1;
    }

    ## Removes ".file_ending" in filename.FILENDING(.gz)
    $file_info_href->{human_genome_reference_name_prefix} =
      fileparse( $human_genome_reference, qr/[.]fasta | [.]fasta[.]gz/xsm );

    $file_info_href->{human_genome_compressed} = $constraint{is_gzipped}->($human_genome_reference);

    if ( $file_info_href->{human_genome_compressed} ) {

        ## Set build file to true to allow for uncompression before analysis
        set_parameter_build_file_status(
            {
                parameter_href => $parameter_href,
                parameter_name => q{human_genome_reference_file_endings},
                status         => 1,
            }
        );
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

sub set_is_sample_files_compressed {

## Function : Set sample files compression status
## Returns  :
## Arguments: $compression_status => Compression status to set
##          : $file_info_href     => File info hash {REF}
##          : $sample_id          => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $compression_status;
    my $file_info_href;
    my $sample_id;

    my $tmpl = {
        compression_status => {
            allow       => [ 0, 1 ],
            defined     => 1,
            required    => 1,
            store       => \$compression_status,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
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

    $file_info_href->{is_files_compressed}{$sample_id} = $compression_status;
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
            allow       => [qw{ in out temp }],
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

        my ( $file_name_prefix, $dirs, $suffix ) = fileparse( $file_path, qr/([.][^.]*)*/sxm );

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
    my ( $infile_path_volume, $file_path_directory ) =
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

    return if ( not $temp_directory );

    ## Also set the temporary file features for stream
    ## Switch to temp dir for path
    my @file_paths_temp =
      map { catfile( $temp_directory, $_ ) } @{ $io_recipe_href->{$stream}{file_names} };

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

sub set_sample_file_attribute {

## Function : Set sample file attributes
## Returns  :
## Arguments: $attribute       => Attribute key
##          : $attribute_value => Attribute value
##          : $file_info_href  => File info hash {REF}
##          : $file_name       => File name
##          : $sample_id       => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $attribute;
    my $attribute_value;
    my $file_info_href;
    my $file_name;
    my $sample_id;

    my $tmpl = {
        attribute => {
            defined     => 1,
            required    => 1,
            store       => \$attribute,
            strict_type => 1,
        },
        attribute_value => {
            required    => 1,
            store       => \$attribute_value,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        file_name => {
            defined     => 1,
            required    => 1,
            store       => \$file_name,
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

    ## Return if nothing to set
    return if ( not defined $attribute_value );

    $file_info_href->{$sample_id}{$file_name}{$attribute} = $attribute_value;
    return;
}

sub set_sample_max_parallel_processes_count {

## Function : Set sample max parallel processes count
## Returns  :
## Arguments: $file_info_href               => File info hash {REF}
##          : $max_parallel_processes_count => New parallel processes count
##          : $sample_id                    => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $max_parallel_processes_count;
    my $sample_id;

    my $tmpl = {
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        max_parallel_processes_count => {
            allow       => qr{ \A \d+ \z }sxm,
            defined     => 1,
            required    => 1,
            store       => \$max_parallel_processes_count,
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

    $file_info_href->{max_parallel_processes_count}{$sample_id} = $max_parallel_processes_count;
    return;
}

sub set_select_file_contigs {

## Function : Set select file contigs
## Returns  :
## Arguments: $file_info_href          => File info hash {REF}
##          : $select_file_contigs_ref => Primary contig hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $select_file_contigs_ref;

    my $tmpl = {
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        select_file_contigs_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$select_file_contigs_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Set select contig sets
    @{ $file_info_href->{select_file_contigs} } = @{$select_file_contigs_ref};

    return;
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
            allow       => [qw{ in out temp }],
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
      @{ $file_info_href->{io}{$upstream_chain_id}{$id}{$upstream_recipe}{$upstream_direction}
          {file_paths} };

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

sub _set_io_files_constant {

## Function : Set the io files per chain and stream for constant features
## Returns  :
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
            allow       => [qw{ in out temp }],
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

  FILE_FEATURE:
    while ( my ( $file_feature, $file_constant_feature ) = each %constant_map ) {

        ## Get unique suffixes
        my @uniq_elements =
          uniq( @{ $io_recipe_href->{$stream}{$file_feature} } );

        next FILE_FEATURE if ( not scalar @uniq_elements == 1 );

        ## Unique - Set file constant suffix
        $io_recipe_href->{$stream}{$file_constant_feature} =
          $uniq_elements[0];
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
            allow       => [qw{ in out temp }],
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
        while ( my ( $file_index, $suffix ) = each @{ $io_recipe_href->{$stream}{file_suffixes} } )
        {

            my $file = $io_recipe_href->{$stream}{$array_key}[$file_index];

            ## Find iterator string between dots
            my ($iterator) = $suffix =~ /[.]([^.]+)[.]/sxm;

            next SUFFIX if ( not $iterator );

            $io_recipe_href->{$stream}{$hash_key}{$iterator} = $file;
        }
    }
    return;
}

1;
