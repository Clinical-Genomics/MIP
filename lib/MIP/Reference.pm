package MIP::Reference;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ fileparse };
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use List::MoreUtils qw { all any };

## MIPs lib/
use MIP::Constants qw{ $COLON $COMMA $LOG_NAME $NEWLINE $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.07;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      check_exome_target_bed_suffix
      check_human_genome_file_endings
      check_nist_file_name
      check_nist_nist_id
      check_nist_sample_id
      check_nist_version
      check_toml_config_for_vcf_tags
      get_dict_contigs
      get_nist_file
      get_select_file_contigs
      parse_meta_file_suffixes
      parse_nist_parameters
      parse_nist_files
      parse_exome_target_bed
      set_nist_file_name_path
      write_contigs_size_file
    };
}

sub check_exome_target_bed_suffix {

## Function : Check that supplied exome target file ends with ".bed" or exit
## Returns  :
## Arguments: $path => Path to check for ".bed" file ending

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $path;

    my $tmpl =
      { path => { defined => 1, required => 1, store => \$path, strict_type => 1, }, };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    if ( $path !~ m{[.]bed$}xsm ) {

        $log->fatal(
            q{Could not find intendended '.bed file ending' for target file: }
              . $path
              . q{ in parameter '--exome_target_bed'},
            $NEWLINE
        );
        exit 1;
    }
    return 1;
}

sub check_human_genome_file_endings {

## Function : Check the existance of associated human genome files
## Returns  :
## Arguments: $human_genome_reference_file_endings_ref => Human genome reference file endings
##          : $human_genome_reference_path             => Human genome reference file path
##          : $parameter_href                          => Parameter hash {REF}
##          : $parameter_name                          => The parameter under evaluation

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $human_genome_reference_file_endings_ref;
    my $human_genome_reference_path;
    my $parameter_href;
    my $parameter_name;

    my $tmpl = {
        human_genome_reference_file_endings_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$human_genome_reference_file_endings_ref,
            strict_type => 1,
        },
        human_genome_reference_path => {
            defined     => 1,
            required    => 1,
            store       => \$human_genome_reference_path,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        parameter_name => { store => \$parameter_name, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Path qw{ check_filesystem_objects_existance };
    use MIP::Parameter qw{ set_parameter_build_file_status };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Count the number of files that exists
    my $existence_check_counter = 0;

  FILE_ENDING:
    foreach my $file_ending ( @{$human_genome_reference_file_endings_ref} ) {

        ## Duplicate to keep fresh with each loop
        my $path = $human_genome_reference_path;

        ## Dict requires no fastq(.gz) ending
        if ( $file_ending eq q{.dict} ) {

            ## Removes ".file_ending" in filename.FILENDING(.gz)
            my ( $file, $dir_path ) =
              fileparse( $path, qr/ [.]fasta | [.]fasta[.]gz /sxm );
            $path = catfile( $dir_path, $file );
        }

        ## Add current ending
        $path = $path . $file_ending;

        my ($does_exist) = check_filesystem_objects_existance(
            {
                object_name    => $path,
                object_type    => q{file},
                parameter_name => $parameter_name,
            }
        );

        ## Sum up the number of file that exists
        $existence_check_counter = $existence_check_counter + $does_exist;
    }
    ## Files need to be built
    if ( $existence_check_counter != scalar @{$human_genome_reference_file_endings_ref} )
    {

        set_parameter_build_file_status {
            (
                parameter_href => $parameter_href,
                parameter_name => $parameter_name,
                status         => 1,
            )
        };
        return;
    }

    # All files exist in this check
    set_parameter_build_file_status {
        (
            parameter_href => $parameter_href,
            parameter_name => $parameter_name,
            status         => 0,
        )
    };
    return;
}

sub check_nist_file_name {

## Function : Check nist file name is defined in nist parameters
## Returns  : 1
## Arguments: $file_name      => File name to check
##          : $nist_id        => Nist id
##          : $nist_parameter => Nist parameter to check
##          : $nist_version   => Nist version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_name;
    my $nist_id;
    my $nist_parameter;
    my $nist_version;

    my $tmpl = {
        file_name => {
            required    => 1,
            store       => \$file_name,
            strict_type => 1,
        },
        nist_id => {
            defined     => 1,
            required    => 1,
            store       => \$nist_id,
            strict_type => 1,
        },
        nist_parameter => {
            defined     => 1,
            required    => 1,
            store       => \$nist_parameter,
            strict_type => 1,
        },
        nist_version => {
            defined     => 1,
            required    => 1,
            store       => \$nist_version,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Require that a file name is defined
    if ( not defined $file_name ) {

        $log->fatal(
            q{Please supply a file name for option: } . join q{=>},
            ( $nist_parameter, $nist_version, $nist_id )
        );
        exit 1;
    }
    return 1;
}

sub check_nist_nist_id {

## Function : Check nist_ids hash contains supplied nist_id in nist parameters
## Returns  : 1
## Arguments: $active_parameter_href => Holds all set parameter for analysis
##          : $nist_id_href          => Nist ids
##          : $nist_parameters_ref   => Nist parameters to check

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $nist_id_href;
    my $nist_parameters_ref;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        nist_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$nist_id_href,
            strict_type => 1,
        },
        nist_parameters_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$nist_parameters_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my %seen;

  NIST_ID:
    foreach my $nist_id ( values %{$nist_id_href} ) {

      NIST_PARAMETER:
        foreach my $nist_parameter ( @{$nist_parameters_ref} ) {

            # Alias
            my $nist_href = \%{ $active_parameter_href->{$nist_parameter} };

          NIST_VERSION:
            foreach my $nist_version ( keys %{$nist_href} ) {

                next NIST_VERSION
                  if ( not exists $nist_href->{$nist_version}{$nist_id} );
                $seen{$nist_id}++;
            }
        }
        next NIST_ID if ( $seen{$nist_id} );

        $log->fatal(
q{Supplied nist id:  $nist_id for option --nist_id is not a defined nist_id supplied in nist options: }
              . join $SPACE,
            @{$nist_parameters_ref}
        );
        exit 1;
    }
    return 1;
}

sub check_nist_sample_id {

## Function : Check nist_ids contain supplied sample_ids
## Returns  : 1
## Arguments: $nist_id_href   => Nist ids
##          : $sample_ids_ref => Sample ids

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $nist_id_href;
    my $sample_ids_ref;

    my $tmpl = {
        nist_id_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$nist_id_href,
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

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

  SAMPLE_ID:
    foreach my $sample_id ( keys %{$nist_id_href} ) {

        ## Find supplied sample id
        next SAMPLE_ID if ( any { $_ eq $sample_id } @{$sample_ids_ref} );

        $log->fatal(
            q{Supplied sample id for option --nist_id }
              . q{ is not one of the included sample ids: }
              . join $SPACE,
            @{$sample_ids_ref}
        );
        exit 1;
    }
    return 1;
}

sub check_nist_version {

## Function : Check nist_versions contain supplied nist_version in nist parameters
## Returns  : 1
## Arguments: $active_parameter_href => Holds all set parameter for analysis
##          : $nist_parameters_ref   => Nist parameters to check

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $nist_parameters_ref;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        nist_parameters_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$nist_parameters_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack
    my @nist_versions = @{ $active_parameter_href->{nist_versions} };

  NIST_PARAMETER:
    foreach my $nist_parameter ( @{$nist_parameters_ref} ) {

        # Alias
        my $nist_href = \%{ $active_parameter_href->{$nist_parameter} };

        ## Check that version exists in nist hashes
        next NIST_PARAMETER
          if ( all { exists $nist_href->{$_} } @nist_versions );

        $log->fatal(
            q{One or more nist versions }
              . ( join $SPACE, @nist_versions )
              . q{ does not exist in nist parameters: }
              . join $SPACE,
            @{$nist_parameters_ref}
        );
        exit 1;
    }
    return 1;
}

sub check_toml_config_for_vcf_tags {

## Function : Check that the toml config contains all neccessary annotation tags
## Returns  : %preop_annotations
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Active_parameter qw{ get_binary_path };
    use MIP::Io::Read qw{ read_from_file };
    use MIP::Vcfanno qw{ check_toml_annotation_for_tags };

    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my $bcftools_binary_path = get_binary_path(
        {
            active_parameter_href => $active_parameter_href,
            binary                => q{bcftools},
        }
    );

    ## TOML parameters
    my %toml = (
        sv_annotate        => q{sv_vcfanno_config},
        variant_annotation => q{vcfanno_config},
    );

    my %missing_tag;
    my %preop_annotations;

  VCFANNO_RECIPE:
    while ( my ( $vcfanno_recipe, $vcfanno_toml_config ) = each %toml ) {

        next VCFANNO_RECIPE if ( not $active_parameter_href->{$vcfanno_recipe} );

        my %vcfanno_config = read_from_file(
            {
                format => q{toml},
                path   => $active_parameter_href->{$vcfanno_toml_config},
            }
        );

      ANNOTATION:
        foreach my $annotation_href ( @{ $vcfanno_config{annotation} } ) {

            ## Only check vcf files
            next ANNOTATION if ( not exists $annotation_href->{fields} );

            ## Check if vcf contains the VCF tag specified in the annotation
            check_toml_annotation_for_tags(
                {
                    annotation_href      => $annotation_href,
                    bcftools_binary_path => $bcftools_binary_path,
                    missing_tag_href     => \%missing_tag,
                    preops_href          => \%preop_annotations,
                }
            );
        }
    }

    if (%missing_tag) {
      REFERENCE_FILE:
        foreach my $reference_file ( keys %missing_tag ) {

            $log->fatal(
                $reference_file
                  . $SPACE
                  . q{misses required ID(s)}
                  . $COLON
                  . $SPACE
                  . join $SPACE,
                @{ $missing_tag{$reference_file} }
            );
        }
        exit 1;
    }

    return %preop_annotations;
}

sub get_dict_contigs {

## Function : Collects sequence contigs used in analysis from human genome sequence
##          : dictionnary (.dict file)
## Returns  : @contigs
## Arguments: $dict_file_path => Dict file path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $dict_file_path;

    my $tmpl = {
        dict_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$dict_file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Environment::Child_process qw{ child_process };
    use MIP::Language::Perl qw{ perl_nae_oneliners };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Build regexp to find contig names
    my @perl_commands = perl_nae_oneliners(
        {
            oneliner_name  => q{get_dict_contigs},
            stdinfile_path => $dict_file_path,
        }
    );

    my @get_dict_contigs_cmds = join $SPACE, ( @perl_commands, );

    # System call
    my %return = child_process(
        {
            commands_ref => \@get_dict_contigs_cmds,
            process_type => q{ipc_cmd_run},
        }
    );

    # Save contigs
    my @contigs = split $COMMA, join $COMMA, @{ $return{stdouts_ref} };

    #my @contigs = split $COMMA, join $COMMA, @{$stdout_buf_ref};

    return @contigs if (@contigs);

    $log->fatal(
        q{Could not detect any 'SN:contig_names' in dict file: } . $dict_file_path );
    exit 1;
}

sub get_nist_file {

## Function : Get nist file
## Returns  : 1
## Arguments: $nist_href     => Nist hash to set file path in
##          : $nist_id       => Nist id
##          : $nist_version  => Nist version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $nist_href;
    my $nist_id;
    my $nist_version;

    my $tmpl = {
        nist_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$nist_href,
            strict_type => 1,
        },
        nist_id => {
            defined     => 1,
            required    => 1,
            store       => \$nist_id,
            strict_type => 1,
        },
        nist_version => {
            defined     => 1,
            required    => 1,
            store       => \$nist_version,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Return nist file
    return $nist_href->{$nist_version}{$nist_id};
}

sub get_select_file_contigs {

## Function : Collects sequences contigs used in select file
## Returns  : @contigs
## Arguments: $select_file_path => Select file path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $select_file_path;

    my $tmpl = {
        select_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$select_file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Environment::Child_process qw{ child_process };
    use MIP::Language::Perl qw{ perl_nae_oneliners };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Build regexp to get contig names
    my @get_select_contigs_cmds = perl_nae_oneliners(
        {
            oneliner_name  => q{get_select_contigs_by_col},
            stdinfile_path => $select_file_path,
        }
    );

    # System call
    my %process_return = child_process(
        {
            commands_ref => \@get_select_contigs_cmds,
            process_type => q{open3},
        }
    );

    # Save contigs
    my @contigs = split $COMMA, join $COMMA, @{ $process_return{stdouts_ref} };

    if ( not @contigs ) {

        $log->fatal(
            q{Could not detect any '##contig' in meta data header in select file: }
              . $select_file_path );
        exit 1;
    }
    return @contigs;
}

sub parse_meta_file_suffixes {

## Function : Checks files to be built by combining object name prefix with suffix.
## Returns  :
## Arguments: $active_parameter_href  => Active parameters for this analysis hash {REF}
##          : $file_name              => File name
##          : $meta_file_suffixes_ref => Reference to the meta file suffixes to be added to the file name {REF}
##          : $parameter_href         => Parameter hash {REF}
##          : $parameter_name         => MIP parameter name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_name;
    my $meta_file_suffixes_ref;
    my $parameter_href;
    my $parameter_name;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        file_name => {
            defined     => 1,
            required    => 1,
            store       => \$file_name,
            strict_type => 1,
        },
        meta_file_suffixes_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$meta_file_suffixes_ref,
            strict_type => 1,
        },
        parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$parameter_href,
            strict_type => 1,
        },
        parameter_name => {
            defined     => 1,
            required    => 1,
            store       => \$parameter_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Path qw{ check_filesystem_objects_existance };
    use MIP::Parameter qw{ get_parameter_attribute set_parameter_build_file_status };

    ## Count the number of files that exists
    my $existence_check_counter = 0;

    my $build_status = 0;

    ## Get parameter object type i.e file or directory
    my $system_object_type = get_parameter_attribute(
        {
            attribute      => q{exists_check},
            parameter_href => $parameter_href,
            parameter_name => $parameter_name,
        }
    );

  FILE_SUFFIX:
    foreach my $file_suffix ( @{$meta_file_suffixes_ref} ) {

        my ($exist) = check_filesystem_objects_existance(
            {
                object_name    => catfile( $file_name . $file_suffix ),
                object_type    => $system_object_type,
                parameter_name => $parameter_name,
            }
        );
        ## Sum up the number of file that exists
        $existence_check_counter = $existence_check_counter + $exist;
    }

    ## Files need to be built
    if ( $existence_check_counter != scalar @{$meta_file_suffixes_ref} ) {

        $build_status = 1;
    }

    # Set build status for parameter
    set_parameter_build_file_status(
        {
            parameter_href => $parameter_href,
            parameter_name => $parameter_name,
            status         => $build_status,
        }
    );
    return;
}

sub parse_nist_files {

## Function : Parse nist files
## Returns  : 1
## Arguments: $active_parameter_href => Holds all set parameter for analysis
##          : $nist_parameters_ref   => Nist parameters to check

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $nist_parameters_ref;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        nist_parameters_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$nist_parameters_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Path qw{ check_filesystem_objects_and_index_existance };

    ## Unpack
    my $reference_dir = $active_parameter_href->{reference_dir};

  NIST_PARAMETER:
    foreach my $nist_parameter ( @{$nist_parameters_ref} ) {

        # Alias
        my $nist_href = \%{ $active_parameter_href->{$nist_parameter} };

      NIST_VERSION:
        foreach my $nist_version ( keys %{$nist_href} ) {

          NIST_FILE:
            while ( my ( $nist_id, $file_name ) = each %{ $nist_href->{$nist_version} } )
            {

                check_nist_file_name(
                    {
                        file_name      => $file_name,
                        nist_id        => $nist_id,
                        nist_parameter => $nist_parameter,
                        nist_version   => $nist_version,
                    }
                );

                set_nist_file_name_path(
                    {
                        file_name     => $file_name,
                        nist_href     => $nist_href,
                        nist_id       => $nist_id,
                        nist_version  => $nist_version,
                        reference_dir => $reference_dir,
                    }
                );

                my $nist_file_path = get_nist_file(
                    {
                        nist_href    => $nist_href,
                        nist_id      => $nist_id,
                        nist_version => $nist_version,
                    }
                );
                ## Check path object exists
                check_filesystem_objects_and_index_existance(
                    {
                        object_name    => ( join q{=>}, ( $nist_version, $nist_id ) ),
                        object_type    => q{file},
                        parameter_name => $nist_parameter,
                        path           => $nist_file_path,
                    }
                );
            }
        }
    }
    return 1;
}

sub parse_nist_parameters {

## Function : Parse nist parameters. Check and add reference directory to file_names.
## Returns  : 1
## Arguments: $active_parameter_href => Holds all set parameter for analysis

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return
      if ( not $active_parameter_href->{rtg_vcfeval} );

    # Unpack
    my %nist_id         = %{ $active_parameter_href->{nist_id} };
    my @nist_parameters = (qw{ nist_call_set_vcf nist_call_set_bed });

    ## Check nist_ids contain supplied sample_ids
    check_nist_sample_id(
        {
            nist_id_href   => \%nist_id,
            sample_ids_ref => $active_parameter_href->{sample_ids},
        }
    );

    ## Check nist_ids contain supplied nist_id in nist parameter
    check_nist_nist_id(
        {
            active_parameter_href => $active_parameter_href,
            nist_id_href          => \%nist_id,
            nist_parameters_ref   => \@nist_parameters,
        }
    );

    ## Check nist_versions contain supplied nist_version in nist parameters
    check_nist_version(
        {
            active_parameter_href => $active_parameter_href,
            nist_parameters_ref   => \@nist_parameters,
        }
    );

    ## Parse nist files
    parse_nist_files(
        {
            active_parameter_href => $active_parameter_href,
            nist_parameters_ref   => \@nist_parameters,
        }
    );

    return 1;
}

sub parse_exome_target_bed {

## Function : Update exome_target_bed files with human genome reference source and version.
##          : Check for correct file suffix
## Returns  :
## Arguments: $exome_target_bed_file_href     => Exome target bed {REF}
##          : $human_genome_reference_source  => Human genome reference source
##          : $human_genome_reference_version => Human genome reference version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $exome_target_bed_file_href;
    my $human_genome_reference_source;
    my $human_genome_reference_version;

    my $tmpl = {
        exome_target_bed_file_href =>
          { required => 1, store => \$exome_target_bed_file_href, },
        human_genome_reference_source => {
            defined     => 1,
            required    => 1,
            store       => \$human_genome_reference_source,
            strict_type => 1,
        },
        human_genome_reference_version => {
            defined     => 1,
            required    => 1,
            store       => \$human_genome_reference_version,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  EXOME_FILE:
    foreach my $exome_target_bed_file ( keys %{$exome_target_bed_file_href} ) {

        my $original_file_name = $exome_target_bed_file;

        ## Replace with actual version
        if ( $exome_target_bed_file =~
            s/genome_reference_source/$human_genome_reference_source/xsm
            && $exome_target_bed_file =~ s/_version/$human_genome_reference_version/xsm )
        {

            ## The delete operator returns the value being deleted
            ## i.e. updating hash key while preserving original info
            $exome_target_bed_file_href->{$exome_target_bed_file} =
              delete $exome_target_bed_file_href->{$original_file_name};
        }

        ## Check that supplied target file ends with ".bed" and otherwise croaks
        check_exome_target_bed_suffix(
            {
                path => $exome_target_bed_file,
            }
        );
    }
    return;
}

sub set_nist_file_name_path {

## Function : Set nist file name path by adding reference directory
## Returns  : 1
## Arguments: $file_name     => File name to prepend reference dir to
##          : $nist_href     => Nist hash to set file path in
##          : $nist_id       => Nist id
##          : $nist_version  => Nist version
##          : $reference_dir => Reference dir path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_name;
    my $nist_href;
    my $nist_id;
    my $nist_version;
    my $reference_dir;

    my $tmpl = {
        file_name => {
            defined     => 1,
            required    => 1,
            store       => \$file_name,
            strict_type => 1,
        },
        nist_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$nist_href,
            strict_type => 1,
        },
        nist_id => {
            defined     => 1,
            required    => 1,
            store       => \$nist_id,
            strict_type => 1,
        },
        nist_version => {
            defined     => 1,
            required    => 1,
            store       => \$nist_version,
            strict_type => 1,
        },
        reference_dir => => {
            defined     => 1,
            required    => 1,
            store       => \$reference_dir,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Add reference directory to path
    $nist_href->{$nist_version}{$nist_id} =
      catfile( $reference_dir, $file_name );
    return 1;
}

sub write_contigs_size_file {

## Function : Write contig size file from human genome sequence fai (.fai) file
## Returns  :
## Arguments: $fai_file_path => Fai file path
##          : $outfile_path  => Chromosome size file path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $fai_file_path;
    my $outfile_path;

    my $tmpl = {
        fai_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$fai_file_path,
            strict_type => 1,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Environment::Child_process qw{ child_process };
    use MIP::Language::Perl qw{ perl_nae_oneliners };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Build regexp to find contig names
    my @perl_commands = perl_nae_oneliners(
        {
            oneliner_name   => q{write_contigs_size_file},
            stdinfile_path  => $fai_file_path,
            stdoutfile_path => $outfile_path,
        }
    );

    my @write_contigs_size_cmd = join $SPACE, ( @perl_commands, );

    # System call
    my %process_return = child_process(
        {
            commands_ref => \@write_contigs_size_cmd,
            process_type => q{ipc_cmd_run},
            verbose      => 0,
        }
    );

    return if ( not @{ $process_return{stderrs_ref} } );

    $log->fatal(q{Could not write contigs size file});
    $log->fatal( q{Error: } . join $NEWLINE, @{ $process_return{stderrs_ref} } );
    exit 1;
}

1;
