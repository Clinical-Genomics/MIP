package MIP::Sample_info;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use Time::Piece;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

## MIPs lib/
use MIP::Constants qw{ $COLON $DOT $EMPTY_STR $LOG_NAME $NEWLINE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{Exporter};

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      get_case_members_attributes_in_duos
      get_family_member_id
      get_read_group
      get_rg_header_line
      get_pedigree_sample_id_attributes
      get_sample_info_case_recipe_attributes
      get_sample_info_sample_recipe_attributes
      reload_previous_pedigree_info
      set_file_path_to_store
      set_gene_panel
      set_in_sample_info
      set_infile_info
      set_no_dry_run_parameters
      set_no_read_direction_file_attributes
      set_parameter_in_sample_info
      set_processing_metafile_in_sample_info
      set_read_direction_file_attributes
      set_recipe_metafile_in_sample_info
      set_recipe_outfile_in_sample_info
      set_sample_gender
      write_sample_info_to_file
    };

}

sub get_case_members_attributes_in_duos {

## Function : Get the sample IDs and phenotypes of the family members
## Returns  : %case_members_attributes
## Arguments: $sample_info_href => Sample info hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_info_href;

    my $tmpl = {
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %case_members_attributes = (
        father   => 0,
        mother   => 0,
        children => [],
        affected => [],
        unknown  => [],
    );

  SAMPLE_ID:
    foreach my $sample_id ( keys %{ $sample_info_href->{sample} } ) {

        my %sample_attributes = get_pedigree_sample_id_attributes(
            {
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        my $phenotype = $sample_attributes{phenotype};
        push @{ $case_members_attributes{$phenotype} }, $sample_id;

        next SAMPLE_ID if ( not( $sample_attributes{father} or $sample_attributes{mother} ) );

        ## Append child
        push @{ $case_members_attributes{children} }, $sample_id;

        $case_members_attributes{father} = $sample_attributes{father};
        $case_members_attributes{mother} = $sample_attributes{mother};
    }

    return %case_members_attributes;
}

sub get_family_member_id {

## Function : Get the sample IDs of the family members
## Returns  : %family_member_id
## Arguments: $sample_info_href => Sample info hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_info_href;

    my $tmpl = {
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %family_member_id = (
        father            => 0,
        mother            => 0,
        children          => [],
        affected          => [],
        unknown_phenotype => [],
    );

  SAMPLE_ID:
    foreach my $sample_id ( keys %{ $sample_info_href->{sample} } ) {

        my $mother = get_pedigree_sample_id_attributes(
            {
                attribute        => q{mother},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        my $father = get_pedigree_sample_id_attributes(
            {
                attribute        => q{father},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        my $phenotype = get_pedigree_sample_id_attributes(
            {
                attribute        => q{phenotype},
                sample_id        => $sample_id,
                sample_info_href => $sample_info_href,
            }
        );

        if ( $phenotype eq q{affected} ) {

            push @{ $family_member_id{affected} }, $sample_id;
        }
        elsif ( $phenotype eq q{unknown} ) {

            push @{ $family_member_id{unknown_phenotype} }, $sample_id;
        }

        next if ( not $father or not $mother );

        ## Append child
        push @{ $family_member_id{children} }, $sample_id;

        $family_member_id{father} = $father;
        $family_member_id{mother} = $mother;
    }

    return %family_member_id;
}

sub get_pedigree_sample_id_attributes {

## Function : Get pedigree sample id attribute
## Returns  : $attribute | %attribute
## Arguments: $attribute        => Attribute key
##          : $sample_id        => Sample id to get attribute for
##          : $sample_info_href => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $attribute;
    my $sample_id;
    my $sample_info_href;

    my $tmpl = {
        attribute => {
            allow => [
                qw{ analysis_type
                  capture_kit
                  dna_sample_id
                  expected_coverage
                  father
                  mother
                  phenotype
                  sample_display_name
                  sample_id
                  sex
                  time_point
                  }
            ],
            store       => \$attribute,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Data::Diver qw{ Dive };

    if ( not $attribute ) {

        ## Return entire sample id hash
        return %{ $sample_info_href->{sample}{$sample_id} };
    }
    ## Get attribute
    if ( defined Dive( $sample_info_href, ( q{sample}, $sample_id, $attribute ) ) ) {

        return $sample_info_href->{sample}{$sample_id}{$attribute};
    }
    return;
}

sub get_read_group {

## Function : Builds hash with read group headers
## Returns  : %read_group
## Arguments: $infile_prefix    => Name of Fastq file minus read direction information
##          : $platform         => Sequencing platform
##          : $sample_id        => Sample ID
##          : $sample_info_href => Sample info hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_prefix;
    my $platform;
    my $sample_id;
    my $sample_info_href;

    my $tmpl = {
        infile_prefix => {
            defined     => 1,
            required    => 1,
            store       => \$infile_prefix,
            strict_type => 1,
        },
        platform => {
            defined     => 1,
            required    => 1,
            store       => \$platform,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Traverse down $sample_info_href
    my %fastq_file =
      %{ $sample_info_href->{sample}{$sample_id}{file}{$infile_prefix}
          {read_direction_file}{ $infile_prefix . q{_1} } };

    my $platform_unit =
      $fastq_file{flowcell} . $DOT . $fastq_file{lane} . $DOT . $fastq_file{sample_barcode};

    ## RG hash
    my %rg = (
        id   => $infile_prefix,
        lane => $fastq_file{lane},
        lb => $sample_id, ## Add molecular library (Dummy value since the actual LB isn't available)
        pl => $platform,
        pu => $platform_unit,
        sm => $sample_id,
    );

    return %rg;
}

sub get_rg_header_line {

## Function : Builds line using read group headers
## Returns  : $rg_header_line
## Arguments: $infile_prefix    => Name of fastq file minus read direction information
##          : $platform         => Sequencing platform
##          : $sample_id        => Sample ID
##          : $sample_info_href => Sample info hash {REF}
##          : $separator        => Separator of the read group elements

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_prefix;
    my $platform;
    my $sample_id;
    my $sample_info_href;
    my $separator;

    my $tmpl = {
        infile_prefix => {
            defined     => 1,
            required    => 1,
            store       => \$infile_prefix,
            strict_type => 1,
        },
        platform => {
            defined     => 1,
            required    => 1,
            store       => \$platform,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        separator => {
            defined     => 1,
            required    => 1,
            store       => \$separator,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %rg = get_read_group(
        {
            infile_prefix    => $infile_prefix,
            platform         => $platform,
            sample_id        => $sample_id,
            sample_info_href => $sample_info_href,
        }
    );

    ## Construct read group line;
    my @rg_elements    = map { uc . $COLON . $rg{$_} } qw{ id lb pl pu sm };
    my $rg_header_line = join $separator, @rg_elements;

    return $rg_header_line;
}

sub get_sample_info_case_recipe_attributes {

## Function : Get case recipe attributes from sample_info hash
## Returns  : "$attribute" or "$attribute_href"
## Arguments: $attribute        => Attribute key
##          : $recipe_name      => Recipe to get attributes from
##          : $sample_info_href => Sample info hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $attribute;
    my $recipe_name;
    my $sample_info_href;

    my $tmpl = {
        attribute => {
            allow       => [qw{ outdirectory outfile path version }],
            store       => \$attribute,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Get and return attribute value
    if ( defined $attribute && $attribute ) {

        return $sample_info_href->{recipe}{$recipe_name}{$attribute};
    }

    ## Get recipe attribute hash
    return %{ $sample_info_href->{recipe}{$recipe_name} };
}

sub get_sample_info_sample_recipe_attributes {

## Function : Get sample recipe attributes from sample_info hash
## Returns  : "$attribute" or "$attribute_href"
## Arguments: $attribute        => Attribute key
##          : $infile           => Infile key
##          : $recipe_name      => Recipe to get attributes from
##          : $sample_id        => Sample id
##          : $sample_info_href => Sample info hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $attribute;
    my $infile;
    my $recipe_name;
    my $sample_id;
    my $sample_info_href;

    my $tmpl = {
        attribute => {
            store       => \$attribute,
            strict_type => 1,
        },
        infile => {
            defined     => 1,
            required    => 1,
            store       => \$infile,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Get and return attribute value
    if ( defined $attribute && $attribute ) {

        return $sample_info_href->{sample}{$sample_id}{recipe}{$recipe_name}{$infile}{$attribute};
    }
    if ( ref $sample_info_href->{sample}{$sample_id}{recipe}{$recipe_name}{$infile} eq q{HASH} ) {

## Get recipe attribute for infile hash
        return %{ $sample_info_href->{sample}{$sample_id}{recipe}{$recipe_name}{$infile} };
    }

    ## No infile level
    ## Get recipe attribute for recipe hash
    return %{ $sample_info_href->{sample}{$sample_id}{recipe}{$recipe_name} };
}

sub reload_previous_pedigree_info {

## Function : Updates sample_info hash with previous run pedigree info
## Returns  :
## Arguments: $sample_info_file_path => Previuos sample info file
##          : $sample_info_href      => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_info_file_path;
    my $sample_info_href;

    my $tmpl = {
        sample_info_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$sample_info_file_path,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Io::Read qw{ read_from_file };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    return if ( not -f $sample_info_file_path );

    # Load parameters from sample_info_file from previous run
    my %previous_sample_info = read_from_file(
        {
            format => q{yaml},
            path   => $sample_info_file_path,
        }
    );

    $log->info( q{Loaded: } . $sample_info_file_path, $NEWLINE );

    ## Update sample_info with pedigree information from previous run
    %{$sample_info_href} = _update_sample_info_hash_pedigree_data(
        {
            sample_info_href          => $sample_info_href,
            previous_sample_info_href => \%previous_sample_info,
        }
    );
    return;
}

sub set_file_path_to_store {

## Function : Set file to store in sample_info
## Returns  :
## Arguments: $format           => File format type
##          : $id               => Id associated with file (sample_id|case_id)
##          : $path             => Path of file
##          : $path_index       => Path of file index
##          : $recipe_name      => Recipe name that produced the file
##          : $sample_info_href => Info on samples and case hash {REF}
##          : $tag              => File tag

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $format;
    my $id;
    my $path;
    my $path_index;
    my $recipe_name;
    my $sample_info_href;
    my $tag;

    my $tmpl = {
        format => {
            allow       => [qw{ bam bb bcf bed bw cram fastq meta png tar vcf wig }],
            defined     => 1,
            required    => 1,
            store       => \$format,
            strict_type => 1,
        },
        id => {
            defined     => 1,
            required    => 1,
            store       => \$id,
            strict_type => 1,
        },
        path => {
            defined     => 1,
            required    => 1,
            store       => \$path,
            strict_type => 1,
        },
        path_index => {
            store       => \$path_index,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        tag => {
            store       => \$tag,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Build file meta data hash
    my %file_info = (
        format     => $format,
        id         => $id,
        path       => $path,
        path_index => $path_index,
        step       => $recipe_name,
        tag        => $tag,
    );

    ## Remove old entries with the same path
    @{ $sample_info_href->{files} } =
      grep { $_->{path} ne $path } @{ $sample_info_href->{files} };

    ## Set file path according to file type and tag
    push @{ $sample_info_href->{files} }, {%file_info};

    return;
}

sub set_gene_panel {

## Function : Collect databases(s) from a database file and sets them to sample_info
## Returns  :
## Arguments: $aggregate_gene_panel_file => Database file
##          : $aggregate_gene_panels_key => Database key i.e. select or range
##          : $recipe_name               => Recipe name
##          : $sample_info_href          => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $aggregate_gene_panel_file;
    my $aggregate_gene_panels_key;
    my $recipe_name;
    my $sample_info_href;

    my $tmpl = {
        aggregate_gene_panel_file => { store => \$aggregate_gene_panel_file, strict_type => 1, },
        aggregate_gene_panels_key => {
            defined     => 1,
            required    => 1,
            store       => \$aggregate_gene_panels_key,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };
    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Language::Perl qw{ perl_nae_oneliners };
    use MIP::Environment::Child_process qw{ child_process };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    return if ( not $aggregate_gene_panel_file );

    # Collect each gene panel features
    my @get_gene_panel_header_cmds = perl_nae_oneliners(
        {
            oneliner_name  => q{get_gene_panel_header},
            stdinfile_path => $aggregate_gene_panel_file,
        }
    );

    my %return = child_process(
        {
            commands_ref => \@get_gene_panel_header_cmds,
            process_type => q{open3},
        }
    );

    my $ret = $return{stdouts_ref}[0]
      or $log->logdie(qq{Unable to parse gene panel information from $aggregate_gene_panel_file});

  LINE:

    # Loop over each gene panel meta data header line
    foreach my $line ( split /:/sxm, $ret ) {

        # Split each info line into gene panel hash
        my %gene_panel = ( split /[,=]/sxm, $line );

        # Gene panel name must exist
        if (    ( defined $gene_panel{gene_panel} )
            and ( not $gene_panel{gene_panel} eq $EMPTY_STR ) )
        {

            # Create unique gene panel ID
            my $gene_panel_name = $gene_panel{gene_panel};
            $sample_info_href->{$recipe_name}{$aggregate_gene_panels_key}{gene_panel}
              {$gene_panel_name} = \%gene_panel;
        }
        else {

            $log->warn(
qq{Unable to write $aggregate_gene_panels_key aggregate gene panel(s) to qc_sample_info. Lacking ##gene_panel=<ID=[?] in aggregate gene panel(s) header.}
            );
        }
    }

    # Call set_processing_metafile_in_sample_info with case parameter
    set_processing_metafile_in_sample_info(
        {
            metafile_tag     => $aggregate_gene_panels_key,
            sample_info_href => $sample_info_href,
            path             => $aggregate_gene_panel_file,
        }
    );
    return;
}

sub set_infile_info {

## Function : Sets information derived from fastq infile name or header to file_info and sample_info hash
## Returns  : $lane_tracker
## Arguments: $file_info_href   => File info hash {REF}
##          : $file_name        => File name
##          : $lane_tracker     => Counts the number of lanes sequenced {REF}
##          : $sample_id        => Sample id
##          : $sample_info_href => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_info_href;
    my $file_name;
    my $lane_tracker;
    my $sample_id;
    my $sample_info_href;

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
        lane_tracker => {
            defined     => 1,
            required    => 1,
            store       => \$lane_tracker,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Fastq qw{ define_mip_fastq_file_features };
    use MIP::File_info qw{
      add_sample_no_direction_infile_prefixes
      get_sample_file_attribute
      set_sample_file_attribute
      set_sample_max_parallel_processes_count
    };

    my %attribute = get_sample_file_attribute(
        {
            file_info_href => $file_info_href,
            file_name      => $file_name,
            sample_id      => $sample_id,
        }
    );

    my $parsed_date = Time::Piece->strptime( $attribute{date}, q{%y%m%d} );
    $parsed_date = $parsed_date->ymd;

    my ( $mip_file_format, $mip_file_format_with_direction,
        $original_file_name_prefix, $run_barcode )
      = define_mip_fastq_file_features(
        {
            date               => $attribute{date},
            direction          => $attribute{direction},
            flowcell           => $attribute{flowcell},
            index              => $attribute{index},
            lane               => $attribute{lane},
            original_file_name => $file_name,
            sample_id          => $sample_id,
        }
      );

    if ( $attribute{direction} == 1 ) {

        add_sample_no_direction_infile_prefixes(
            {
                file_info_href  => $file_info_href,
                mip_file_format => $mip_file_format,
                sample_id       => $sample_id,
            }
        );
        my $sequence_run_type = $attribute{is_interleaved} ? q{interleaved} : q{single-end};
        set_sample_file_attribute(
            {
                attribute       => q{sequence_run_type},
                attribute_value => $sequence_run_type,
                file_info_href  => $file_info_href,
                file_name       => $mip_file_format,
                sample_id       => $sample_id,
            }
        );

        set_sample_max_parallel_processes_count(
            {
                file_info_href               => $file_info_href,
                max_parallel_processes_count => $lane_tracker,
                sample_id                    => $sample_id,
            }
        );

        my %direction_one_metric = (
            interleaved       => $attribute{is_interleaved},
            sequence_length   => $attribute{read_length},
            sequence_run_type => q{single-end},
        );
        set_no_read_direction_file_attributes(
            {
                file_name              => $mip_file_format,
                no_read_direction_href => \%direction_one_metric,
                sample_id              => $sample_id,
                sample_info_href       => $sample_info_href,
            }
        );
        $lane_tracker++;
    }
    if ( $attribute{direction} == 2 ) {

        set_sample_file_attribute(
            {
                attribute       => q{sequence_run_type},
                attribute_value => q{paired-end},
                file_info_href  => $file_info_href,
                file_name       => $mip_file_format,
                sample_id       => $sample_id,
            }
        );

        my %direction_two_metric = ( sequence_run_type => q{paired-end}, );
        set_no_read_direction_file_attributes(
            {
                file_name              => $mip_file_format,
                no_read_direction_href => \%direction_two_metric,
                sample_id              => $sample_id,
                sample_info_href       => $sample_info_href,
            }
        );
    }

    my %both_directions_metric = (
        date                      => $parsed_date,
        flowcell                  => $attribute{flowcell},
        lane                      => $attribute{lane},
        original_file_name        => $file_name,
        original_file_name_prefix => $original_file_name_prefix,
        read_direction            => $attribute{direction},
        run_barcode               => $run_barcode,
        sample_barcode            => $attribute{index},
    );

    set_read_direction_file_attributes(
        {
            direction_file_name => $mip_file_format_with_direction,
            file_name           => $mip_file_format,
            read_direction_href => \%both_directions_metric,
            sample_id           => $sample_id,
            sample_info_href    => $sample_info_href,
        }
    );
    return $lane_tracker;
}

sub set_no_dry_run_parameters {

## Function : Set parameters for true run i.e. not a dry run
## Returns  :
## Arguments: $analysis_date    => Analysis date
##          : $is_dry_run_all   => Dry run boolean
##          : $mip_version      => MIP version
##          : $sample_info_href => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $analysis_date;
    my $is_dry_run_all;
    my $mip_version;
    my $sample_info_href;

    my $tmpl = {
        analysis_date => {
            defined     => 1,
            required    => 1,
            store       => \$analysis_date,
            strict_type => 1,
        },
        is_dry_run_all => {
            allow       => [ 0, 1, undef ],
            required    => 1,
            store       => \$is_dry_run_all,
            strict_type => 1,
        },
        mip_version => {
            defined     => 1,
            required    => 1,
            store       => \$mip_version,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return if ($is_dry_run_all);

    my %no_dry_run_info = (
        analysisrunstatus => q{not_finished},
        analysis_date     => $analysis_date,
        mip_version       => $mip_version,
    );

  PARAMETER_NAME:
    while ( my ( $parameter_name, $parameter_value ) = each %no_dry_run_info ) {

        set_in_sample_info(
            {
                key              => $parameter_name,
                sample_info_href => $sample_info_href,
                value            => $parameter_value,
            }
        );
    }
    return;
}

sub set_no_read_direction_file_attributes {

##Function : Sets file attributes for no read direction file
##Returns  :
##Arguments: $file_name              => File name
##         : $no_read_direction_href => No read direction keys and values to set
##         : $sample_id              => Sample id
##         : $sample_info_href       => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $file_name;
    my $no_read_direction_href;
    my $sample_id;
    my $sample_info_href;

    my $tmpl = {
        file_name => {
            defined     => 1,
            required    => 1,
            store       => \$file_name,
            strict_type => 1,
        },
        no_read_direction_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$no_read_direction_href,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    while ( my ( $attribute, $attribute_value ) = each %{$no_read_direction_href} ) {

        $sample_info_href->{sample}{$sample_id}{file}{$file_name}{$attribute} =
          $attribute_value;
    }
    return;
}

sub set_in_sample_info {

##Function : Sets key and value in sample info
##Returns  :
##Arguments: $key              => Key to add
##         : $sample_info_href => Info on samples and case hash {REF}
##         : $value            => Value to add

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $key;
    my $sample_info_href;
    my $value;

    my $tmpl = {
        key => {
            defined     => 1,
            required    => 1,
            store       => \$key,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
        value => {
            required => 1,
            store    => \$value,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return if ( not defined $value );

    $sample_info_href->{$key} = $value;
    return;
}

sub set_read_direction_file_attributes {

##Function : Sets read direction file attributes
##Returns  :
##Arguments: $direction_file_name => File name with read direction
##         : $file_name           => File name
##         : $read_direction_href => Read direction keys and values to set
##         : $sample_id           => Sample id
##         : $sample_info_href    => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $direction_file_name;
    my $file_name;
    my $read_direction_href;
    my $sample_id;
    my $sample_info_href;

    my $tmpl = {
        direction_file_name => {
            defined     => 1,
            required    => 1,
            store       => \$direction_file_name,
            strict_type => 1,
        },
        file_name => {
            defined     => 1,
            required    => 1,
            store       => \$file_name,
            strict_type => 1,
        },
        read_direction_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$read_direction_href,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    while ( my ( $attribute, $attribute_value ) = each %{$read_direction_href} ) {

        $sample_info_href->{sample}{$sample_id}{file}{$file_name}
          {read_direction_file}{$direction_file_name}{$attribute} = $attribute_value;
    }
    return;
}

sub set_recipe_outfile_in_sample_info {

## Function : Sets path and/or outdirectory and/or outfile and/or version from recipes to sample_info to track all outfiles and extract downstream
## Returns  :
## Arguments: $infile           => Infile for data at sample level {Optional}
##          : $outdirectory     => Outdirectory of the file
##          : $outfile          => Outfile name
##          : $path             => Path of file
##          : $recipe_name      => Recipe name
##          : $sample_id        => Sample_id for data at sample level {Optional}
##          : $sample_info_href => Records on samples and case hash {REF}
##          : $version          => Version of file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile;
    my $outdirectory;
    my $outfile;
    my $path;
    my $recipe_name;
    my $sample_id;
    my $sample_info_href;
    my $version;

    my $tmpl = {
        infile       => { strict_type => 1, store => \$infile },
        outdirectory => {
            strict_type => 1,
            store       => \$outdirectory,
        },
        outfile => {
            strict_type => 1,
            store       => \$outfile,
        },
        path => {
            strict_type => 1,
            store       => \$path,
        },
        recipe_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$recipe_name,
        },
        sample_id        => { strict_type => 1, store => \$sample_id, },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        version => {
            strict_type => 1,
            store       => \$version,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Set the key and value pair to add to sample_info hash
    my %parameter = (
        outdirectory => $outdirectory,
        outfile      => $outfile,
        path         => $path,
        version      => $version,
    );

    if ( defined $sample_id && defined $infile ) {

      SAMPLE_PARAMETER:
        while ( my ( $parameter_key, $parameter_value ) = each %parameter ) {

            if ( defined $parameter_value ) {

                $sample_info_href->{sample}{$sample_id}{recipe}{$recipe_name}
                  {$infile}{$parameter_key} = $parameter_value;
            }
        }
    }
    elsif ( defined $sample_id ) {

      SAMPLE_PARAMETER:
        while ( my ( $parameter_key, $parameter_value ) = each %parameter ) {

            if ( defined $parameter_value ) {

                $sample_info_href->{sample}{$sample_id}{recipe}{$recipe_name}{$parameter_key} =
                  $parameter_value;
            }
        }
    }
    else {

      FAMILY_PARAMETER:
        while ( my ( $parameter_key, $parameter_value ) = each %parameter ) {

            if ( defined $parameter_value ) {

                $sample_info_href->{recipe}{$recipe_name}{$parameter_key} =
                  $parameter_value;
            }
        }
    }
    return;
}

sub set_processing_metafile_in_sample_info {

## Function : Sets metafile path from sample_id|case_id processing to sample_info to track all metafiles and extract downstream
## Returns  :
## Arguments: $metafile_tag     => Id tag of meta file
##          : $path             => Path of file
##          : $sample_id        => Sample_id for data at sample level {Optional}
##          : $sample_info_href => Records on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $metafile_tag;
    my $path;
    my $sample_id;
    my $sample_info_href;

    my $tmpl = {
        metafile_tag => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$metafile_tag,
        },
        path => {
            strict_type => 1,
            store       => \$path,
        },
        sample_id        => { strict_type => 1, store => \$sample_id },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Set the key and value pair to add to sample_info hash
    my %parameter = ( path => $path, );

    if ( defined $sample_id ) {

      SAMPLE_PARAMETER:
        while ( my ( $parameter_key, $parameter_value ) = each %parameter ) {

            if ( defined $parameter_value ) {

                $sample_info_href->{sample}{$sample_id}{$metafile_tag}{$parameter_key} =
                  $parameter_value;
            }
        }
    }
    else {

      FAMILY_PARAMETER:
        while ( my ( $parameter_key, $parameter_value ) = each %parameter ) {

            if ( defined $parameter_value ) {

                $sample_info_href->{$metafile_tag}{$parameter_key} = $parameter_value;
            }
        }
    }
    return;
}

sub set_recipe_metafile_in_sample_info {

## Function : Sets path and/or directory and/or file and/or version from recipes to sample_info to track all metafiles and extract downstream
## Returns  :
## Arguments: $infile           => Infile for data at sample level {Optional}
##          : $directory        => Directory of the file
##          : $file             => File name
##          : $metafile_tag     => Id tag of meta file
##          : $path             => Path of file
##          : $processed_by     => Processed by
##          : $recipe_name      => Recipe name
##          : $version          => Version of file
##          : $sample_id        => Sample_id for data at sample level {Optional}
##          : $sample_info_href => Records on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $directory;
    my $file;
    my $metafile_tag;
    my $infile;
    my $path;
    my $processed_by;
    my $recipe_name;
    my $sample_id;
    my $sample_info_href;
    my $version;

    my $tmpl = {
        directory => {
            strict_type => 1,
            store       => \$directory,
        },
        file => {
            strict_type => 1,
            store       => \$file,
        },
        infile       => { strict_type => 1, store => \$infile, },
        metafile_tag => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$metafile_tag,
        },
        path => {
            strict_type => 1,
            store       => \$path,
        },
        processed_by => {
            strict_type => 1,
            store       => \$processed_by,
        },
        recipe_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$recipe_name,
        },
        sample_id        => { strict_type => 1, store => \$sample_id },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
        version => {
            strict_type => 1,
            store       => \$version,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Set the key and value pair to add to sample_info hash
    my %parameter = (
        directory    => $directory,
        file         => $file,
        path         => $path,
        processed_by => $processed_by,
        version      => $version,
    );

    if ( defined $sample_id && defined $infile ) {

      SAMPLE_PARAMETER:
        while ( my ( $parameter_key, $parameter_value ) = each %parameter ) {

            if ( defined $parameter_value ) {

                $sample_info_href->{sample}{$sample_id}{recipe}{$recipe_name}
                  {$infile}{$metafile_tag}{$parameter_key} = $parameter_value;
            }
        }
    }
    else {

      FAMILY_PARAMETER:
        while ( my ( $parameter_key, $parameter_value ) = each %parameter ) {

            if ( defined $parameter_value ) {

                $sample_info_href->{recipe}{$recipe_name}{$metafile_tag}{$parameter_key} =
                  $parameter_value;
            }
        }
    }
    return;
}

sub set_parameter_in_sample_info {

## Function : Sets parameter to sample_info from active_parameter and file_info
## Returns  :
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $file_info_href        => File info hash {REF}
##          : $sample_info_href      => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $sample_info_href;

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
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Pedigree qw{ has_duo has_trio };

    my %set_parameter_map = (
        analysis_type => {
            key    => q{analysis_type},
            set_at => $sample_info_href,
            value  => $active_parameter_href->{analysis_type},
        },
        expected_coverage => {
            key    => q{expected_coverage},
            set_at => $sample_info_href,
            value  => $active_parameter_href->{expected_coverage},
        },
        has_duo => {
            key    => q{has_duo},
            set_at => $sample_info_href,
            value  => has_duo(
                {
                    active_parameter_href => $active_parameter_href,
                    sample_info_href      => $sample_info_href,
                }
            ),
        },
        has_trio => {
            key    => q{has_trio},
            set_at => $sample_info_href,
            value  => has_trio(
                {
                    active_parameter_href => $active_parameter_href,
                    sample_info_href      => $sample_info_href,
                }
            ),
        },
        human_genome_build_path => {
            key    => q{path},
            set_at => \%{ $sample_info_href->{human_genome_build} },
            value  => $active_parameter_href->{human_genome_reference},
        },
        human_genome_build_source => {
            key    => q{source},
            set_at => \%{ $sample_info_href->{human_genome_build} },
            value  => $file_info_href->{human_genome_reference_source},
        },
        human_genome_build_version => {
            key    => q{version},
            set_at => \%{ $sample_info_href->{human_genome_build} },
            value  => $file_info_href->{human_genome_reference_version},
        },
        last_log_file_path => {
            key    => q{last_log_file_path},
            set_at => $sample_info_href,
            value  => $active_parameter_href->{log_file},
        },
        log_file_dir => {
            key    => q{log_file_dir},
            set_at => $sample_info_href,
            value  => dirname( dirname( $active_parameter_href->{log_file} ) ),
        },
        pedigree_file_path => {
            key    => q{path},
            set_at => \%{ $sample_info_href->{pedigree_file} },
            value  => $active_parameter_href->{pedigree_file},
        },
    );

  PARMETER_MAP_HREF:
    foreach my $parameter_href ( values %set_parameter_map ) {

        set_in_sample_info(
            {
                key              => $parameter_href->{key},
                sample_info_href => $parameter_href->{set_at},
                value            => $parameter_href->{value},
            }
        );
    }
    return;
}

sub set_sample_gender {

## Function : Set gender in sample info
## Returns  :
## Arguments: $gender           => Gender
##          : $sample_id        => Sample id
##          : $sample_info_href => Records on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $gender;
    my $sample_id;
    my $sample_info_href;

    my $tmpl = {
        gender => {
            allow       => [qw{ female male other unknown }],
            required    => 1,
            store       => \$gender,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    $sample_info_href->{sample}{$sample_id}{sex} = $gender;

    return;
}

sub write_sample_info_to_file {

## Function : Write sample info to file
## Returns  :
## Arguments: $sample_info_file => Sample info file to write to
##          : $sample_info_href => Records on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_info_file;
    my $sample_info_href;

    my $tmpl = {
        sample_info_file => { strict_type => 1, store => \$sample_info_file },
        sample_info_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$sample_info_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Io::Write qw{ write_to_file };

    return if ( not $sample_info_file );

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Writes a YAML hash to file
    write_to_file(
        {
            data_href => $sample_info_href,
            format    => q{yaml},
            path      => $sample_info_file,
        }
    );
    $log->info( q{Wrote: } . $sample_info_file );

    return;
}

sub _update_sample_info_hash_pedigree_data {

## Function : Update sample_info with information from pedigree from previous run.
##          : Required e.g. if only updating single sample analysis chains from trio.
## Returns  : %{$previous_sample_info_href}
## Arguments: $previous_sample_info_href => Allowed parameters from pedigre file hash {REF}
##          : $sample_info_href          => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $previous_sample_info_href;
    my $sample_info_href;

    my $tmpl = {
        previous_sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$previous_sample_info_href,
            strict_type => 1,
        },
        sample_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$sample_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  SAMPLE_ID:
    foreach my $sample_id ( keys %{ $sample_info_href->{sample} } ) {

        ## Alias
        my $sample_href = \%{ $sample_info_href->{sample}{$sample_id} };
        my $previous_sample_href =
          \%{ $previous_sample_info_href->{sample}{$sample_id} };

      PEDIGREE_KEY:
        foreach my $pedigree_key ( keys %{$sample_href} ) {

            ## Previous run information, which should be updated using pedigree from current analysis
            if ( exists $previous_sample_href->{$pedigree_key} ) {

                ## Required to update keys downstream
                my $previous_pedigree_value = delete $sample_href->{$pedigree_key};

                ## Update previous sample info key
                $previous_sample_href->{$pedigree_key} = $previous_pedigree_value;
                next PEDIGREE_KEY;
            }

            ## New sample_id or key
            $previous_sample_href->{$pedigree_key} = $sample_href->{$pedigree_key};
        }
    }
    return %{$previous_sample_info_href};
}

1;
