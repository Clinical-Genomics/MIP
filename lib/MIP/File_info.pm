package MIP::File_info;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname fileparse };
use File::Spec::Functions qw{ catfile };
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
use MIP::Constants qw{ $GENOME_VERSION $LOG_NAME $NEWLINE $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.12;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      add_sample_fastq_file_lanes
      check_parameter_metafiles
      get_is_sample_files_compressed
      get_sample_file_attribute
      parse_file_compression_features
      parse_files_compression_status
      parse_sample_fastq_file_attributes
      parse_select_file_contigs
      set_alt_loci_contigs
      set_bam_contigs
      set_dict_contigs
      set_file_tag
      set_infiles
      set_is_sample_files_compressed
      set_human_genome_reference_features
      set_primary_contigs
      set_sample_file_attribute
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
##          : $lane           => Fast file name
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
                        active_parameter_href => $active_parameter_href,
                        file_name             => $path,
                        meta_file_suffixes_ref =>
                          \@{ $file_info_href->{$parameter_name} },
                        parameter_href => $parameter_href,
                        parameter_name => $parameter_name,
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
                active_parameter_href => $active_parameter_href,
                file_name             => $active_parameter_href->{human_genome_reference},
                meta_file_suffixes_ref => \@{ $file_info_href->{$parameter_name} },
                parameter_href         => $parameter_href,
                parameter_name         => $parameter_name,
            }
        );
    }
    return;
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

        ## Return entire file name hash
        return %{ $file_info_href->{$sample_id}{$file_name} };
    }
    ## Get attribute
    my $stored_attribute =
      $file_info_href->{$sample_id}{$file_name}{$attribute};

    ## Return requested attribute
    return $stored_attribute;
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

    use MIP::File::Path qw{ check_gzipped };

    my %attribute = (
        is_file_compressed => 0,
        read_file_command  => q{cat},
    );

    ## Check if a file is gzipped.
    my $is_gzipped = check_gzipped( { file_name => $file_name, } );

    ## Gzipped
    if ($is_gzipped) {

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

sub set_bam_contigs {

## Function : Set bam contigs
## Returns  :
## Arguments: $bam_contig_set_name => Bam contig set identifier
##          : $file_info_href      => File info hash {REF}
##          : $primary_contigs_ref => Primary contig hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bam_contig_set_name;
    my $file_info_href;
    my $primary_contigs_ref;

    my $tmpl = {
        bam_contig_set_name => {
            defined     => 1,
            required    => 1,
            store       => \$bam_contig_set_name,
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
    @{ $file_info_href->{$bam_contig_set_name} } = @{$primary_contigs_ref};

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
    use MIP::File::Path qw{ check_gzipped };
    use MIP::Parameter qw{ set_parameter_build_file_status };

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

    $file_info_href->{human_genome_compressed} =
      check_gzipped( { file_name => $human_genome_reference, } );

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

sub set_primary_contigs {

## Function : Set primary contigs
## Returns  :
## Arguments: $file_info_href          => File info hash {REF}
##          : $primary_contigs_ref     => Primary contig hash {REF}
##          : $primary_contig_set_name => Primary contig set identifier

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $primary_contigs_ref;
    my $primary_contig_set_name;

    my $tmpl = {
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
        primary_contig_set_name => {
            defined     => 1,
            required    => 1,
            store       => \$primary_contig_set_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Set primary contig sets
    @{ $file_info_href->{$primary_contig_set_name} } = @{$primary_contigs_ref};

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
1;
