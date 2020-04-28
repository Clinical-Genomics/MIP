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
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $COLON $DOT $EMPTY_STR $LOG_NAME $NEWLINE $SPACE $UNDERSCORE };

BEGIN {

    require Exporter;
    use base qw{Exporter};

    # Set the version for version checking
    our $VERSION = 1.24;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      get_family_member_id
      get_read_group
      get_rg_header_line
      get_pedigree_sample_id_attributes
      get_sample_info_case_recipe_attributes
      get_sample_info_sample_recipe_attributes
      get_sequence_run_type
      get_sequence_run_type_is_interleaved
      reload_previous_pedigree_info
      set_file_path_to_store
      set_gene_panel
      set_infile_info
      set_most_complete_vcf
      set_parameter_in_sample_info
      set_processing_metafile_in_sample_info
      set_recipe_metafile_in_sample_info
      set_recipe_outfile_in_sample_info
      set_in_sample_info
    };

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
## Returns  : $attribute
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

    if ( not $attribute ) {

        ## Return entire sample id hash
        return %{ $sample_info_href->{sample}{$sample_id} };
    }
    ## Get attribute
    my $stored_attribute = $sample_info_href->{sample}{$sample_id}{$attribute};

    ## Return requested attribute
    return $stored_attribute;
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
        $fastq_file{flowcell}
      . $DOT
      . $fastq_file{lane}
      . $DOT
      . $fastq_file{sample_barcode};

    ## RG hash
    my %rg = (
        id   => $infile_prefix,
        lane => $fastq_file{lane},
        lb   => $sample_id
        ,    ## Add molecular library (Dummy value since the actual LB isn't available)
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

        return $sample_info_href->{sample}{$sample_id}{recipe}{$recipe_name}
          {$infile}{$attribute};
    }
    if (
        ref $sample_info_href->{sample}{$sample_id}{recipe}{$recipe_name}{$infile} eq
        q{HASH} )
    {

## Get recipe attribute for infile hash
        return %{ $sample_info_href->{sample}{$sample_id}{recipe}{$recipe_name}{$infile}
        };
    }

    ## No infile level
    ## Get recipe attribute for recipe hash
    return %{ $sample_info_href->{sample}{$sample_id}{recipe}{$recipe_name} };
}

sub get_sequence_run_type {

## Function : Return sequence run type
## Returns  : $sequence_run_type | %sequence_run_type
## Arguments: $infile_lane_prefix      => Infile lane prefix
##          : $infile_lane_prefix_href => Infile lane prefix hash {REF}
##          : $sample_id               => Sample id
##          : $sample_info_href        => Sample info hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_lane_prefix;
    my $infile_lane_prefix_href;
    my $sample_id;
    my $sample_info_href;

    my $tmpl = {
        infile_lane_prefix => {
            store       => \$infile_lane_prefix,
            strict_type => 1,
        },
        infile_lane_prefix_href => {
            default     => {},
            store       => \$infile_lane_prefix_href,
            strict_type => 1,
        },
        sample_id => {
            required    => 1,
            defined     => 1,
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

    if ( %{$infile_lane_prefix_href} ) {

        my %sequence_run_type;

      INFILE_LANE_PREFIX:
        foreach my $infile_lane_prefix ( @{ $infile_lane_prefix_href->{$sample_id} } ) {
            $sequence_run_type{$infile_lane_prefix} =
              $sample_info_href->{sample}{$sample_id}{file}{$infile_lane_prefix}
              {sequence_run_type};
        }

        return %sequence_run_type;
    }
    elsif ( defined $infile_lane_prefix ) {

        return $sample_info_href->{sample}{$sample_id}{file}{$infile_lane_prefix}
          {sequence_run_type};

    }

    croak q{Either $infile_lane_prefix_href or $infile_prefix must be provided!};
}

sub get_sequence_run_type_is_interleaved {

## Function : Get sequence run type interleaved info
## Returns  : undef (not interleaved) | 1 (is interleaved)
## Arguments: $infile_lane_prefix => Infile lane prefix
##          : $sample_id          => Sample id
##          : $sample_info_href   => Sample info hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_lane_prefix;
    my $sample_id;
    my $sample_info_href;

    my $tmpl = {
        infile_lane_prefix => {
            required    => 1,
            defined     => 1,
            store       => \$infile_lane_prefix,
            strict_type => 1,
        },
        sample_id => {
            required    => 1,
            defined     => 1,
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

    return
      if (
        ref $sample_info_href->{sample}{$sample_id}{file}{$infile_lane_prefix}
        {sequence_run_type} ne q{HASH} );

    return $sample_info_href->{sample}{$sample_id}{file}{$infile_lane_prefix}
      {sequence_run_type}{interleaved};

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
            allow       => [qw{ fastq bam bcf bed cram meta tar vcf }],
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

    ## Set file path according to file type and tag
    push @{ $sample_info_href->{files} }, {%file_info};

    return;
}

sub set_gene_panel {

## Function : Collect databases(s) from a database file and sets them to sample_info
## Returns  :
## Arguments: $aggregate_gene_panel_file => Database file
##          : $aggregate_gene_panels_key => Database key i.e. select or range
##          : $case_id                   => Case ID
##          : $log                       => Log object
##          : $recipe_name               => Recipe name
##          : $sample_info_href          => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $aggregate_gene_panel_file;
    my $aggregate_gene_panels_key;
    my $case_id;
    my $log;
    my $recipe_name;
    my $sample_info_href;

    my $tmpl = {
        aggregate_gene_panel_file =>
          { store => \$aggregate_gene_panel_file, strict_type => 1, },
        aggregate_gene_panels_key => {
            defined     => 1,
            required    => 1,
            store       => \$aggregate_gene_panels_key,
            strict_type => 1,
        },
        case_id => {
            defined     => 1,
            required    => 1,
            store       => \$case_id,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
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

    if ($aggregate_gene_panel_file) {

        # Collect each gene panel features
        my %gene_panel;
        my %header = (
            display_name => q{display_name},
            gene_panel   => q{gene_panel},
            updated_at   => q{updated_at},
            version      => q{version},
        );

        # Execute perl
        my $sub_database_regexp = q?perl -nae '?;

        # If line starts with gene panel comment
        $sub_database_regexp .= q?if ($_=~/^##gene_panel=/)? . $SPACE;

        # Remove newline char and split fields
        $sub_database_regexp .= q?{chomp($_);my @entries=split(/,/, $_);? . $SPACE;

        # Join fields with comma separator appending ":". Skip rest if it's a comment
        $sub_database_regexp .=
          q?my $entry = join(",", $_); print $entry.":" } if($_=~/^#\w/) {last;}'?;

        # Collect header_lines(s) from select_file header
        my $ret = `$sub_database_regexp $aggregate_gene_panel_file`;

        # Split each gene panel meta data header line into array element
        my @header_lines = split /:/sxm, $ret;

      LINE:
        foreach my $line (@header_lines) {

            # Split each member database line into features
            my @features = split /,/sxm, $line;

          ELEMENT:
            foreach my $feature_element (@features) {

              KEY_VALUE:
                foreach my $gene_panel_header_element ( keys %header ) {

                    # Parse the features using defined header keys
                    if ( $feature_element =~ /$gene_panel_header_element=/sxm ) {

                        my @temps = split /=/sxm, $feature_element;

                        # Value
                        $gene_panel{ $header{$gene_panel_header_element} } =
                          $temps[1];
                        last;
                    }
                }
            }

            if ( defined $gene_panel{gene_panel} ) {

                # Create unique gene panel ID
                my $gene_panel_name = $gene_panel{gene_panel};

                ## Add new entries
              FEATURE:
                foreach my $feature ( keys %gene_panel ) {

                    $sample_info_href->{$recipe_name}
                      {$aggregate_gene_panels_key}{gene_panel}{$gene_panel_name}{$feature}
                      = $gene_panel{$feature};
                }
            }
            else {

                $log->warn( q{Unable to write}
                      . $SPACE
                      . $aggregate_gene_panels_key
                      . $SPACE
                      . q{aggregate gene panel(s) to qc_sample_info. Lacking ##gene_panel=<ID=[?] or version=[?] in aggregate gene panel(s) header.}
                      . $NEWLINE );
            }

            # Reset hash for next line
            %gene_panel = ();
        }

        # Call set_processing_metafile_in_sample_info with case parameter
        set_processing_metafile_in_sample_info(
            {
                metafile_tag     => $aggregate_gene_panels_key,
                sample_info_href => $sample_info_href,
                path             => $aggregate_gene_panel_file,
            }
        );
    }
    return;
}

sub set_infile_info {

## Function : Sets information derived from infile name to sample_info hash. Tracks the number of lanes sequenced and checks unique array elements.
## Returns  : $lane_tracker
## Arguments: $active_parameter_href           => Active parameters for this analysis hash {REF}
##          : $date                            => Flow-cell sequencing date
##          : $direction                       => Sequencing read direction
##          : $file_index                      => Index of file
##          : $file_info_href                  => File info hash {REF}
##          : $flowcell                        => Flow-cell id
##          : $index                           => The DNA library preparation molecular barcode
##          : $infile_both_strands_prefix_href => The infile(s) without the ".ending" and strand info {REF}
##          : $infile_lane_prefix_href         => Infile(s) without the ".ending" {REF}
##          : $is_interleaved                  => Infile is interleaved
##          : $lane                            => Flow-cell lane
##          : $lane_tracker                    => Counts the number of lanes sequenced {REF}
##          : $read_length                     => Sequence read length
##          : $sample_id                       => Sample id
##          : $sample_info_href                => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $date;
    my $direction;
    my $file_index;
    my $file_info_href;
    my $flowcell;
    my $index;
    my $infile_both_strands_prefix_href;
    my $infile_lane_prefix_href;
    my $is_interleaved;
    my $lane;
    my $lane_tracker;
    my $read_length;
    my $sample_id;
    my $sample_info_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        date      => { defined => 1, required => 1, store => \$date, strict_type => 1, },
        direction => {
            allow       => [ 1, 2 ],
            defined     => 1,
            required    => 1,
            store       => \$direction,
            strict_type => 1,
        },
        file_index => {
            allow       => qr{ \A\d+\z }xsm,
            defined     => 1,
            required    => 1,
            store       => \$file_index,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        flowcell => {
            defined     => 1,
            required    => 1,
            store       => \$flowcell,
            strict_type => 1,
        },
        index => { defined => 1, required => 1, store => \$index, strict_type => 1, },
        infile_both_strands_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_both_strands_prefix_href,
            strict_type => 1,
        },
        infile_lane_prefix_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$infile_lane_prefix_href,
            strict_type => 1,
        },
        is_interleaved => {
            required    => 1,
            store       => \$is_interleaved,
            strict_type => 1,
        },
        lane => {
            allow       => qr{ \A\d+\z }xsm,
            defined     => 1,
            required    => 1,
            store       => \$lane,
            strict_type => 1,
        },
        lane_tracker => {
            defined     => 1,
            required    => 1,
            store       => \$lane_tracker,
            strict_type => 1,
        },
        read_length => {
            required    => 1,
            store       => \$read_length,
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

    my $parsed_date = Time::Piece->strptime( $date, q{%y%m%d} );
    $parsed_date = $parsed_date->ymd;

    my ( $mip_file_format, $mip_file_format_with_direction,
        $original_file_name_prefix, $run_barcode )
      = _file_name_formats(
        {
            date               => $date,
            direction          => $direction,
            flowcell           => $flowcell,
            index              => $index,
            lane               => $lane,
            original_file_name => $file_info_href->{$sample_id}{mip_infiles}[$file_index],
            sample_id          => $sample_id,
        }
      );

    ## Read 1
    if ( $direction == 1 ) {

        ## Add lane
        push @{ $file_info_href->{$sample_id}{lanes} }, $lane;

# Save new format (sample_id_date_flow-cell_index_lane) in hash with samplid as keys and inputfiles in array. Note: These files have not been created yet and there is one entry into hash for both strands and the file suffix is removed (.fastq).
        $infile_lane_prefix_href->{$sample_id}[$lane_tracker] = $mip_file_format;

        ## Detect Undetermined in flowcell id
        if ( $flowcell =~ /Undetermined/ixsm ) {

            ## Set Undetermined to true for file
            $file_info_href->{undetermined_in_file_name}{$mip_file_format} = 1;
        }

        my %direction_one_metric = (
            interleaved       => $is_interleaved,
            sequence_length   => $read_length,
            sequence_run_type => q{single-end},
        );

        ## Alias
        my $file_level_href =
          \%{ $sample_info_href->{sample}{$sample_id}{file}{$mip_file_format} };
      INFO:
        while ( my ( $file_key, $file_value ) = each %direction_one_metric ) {

            $file_level_href->{$file_key} = $file_value;
        }

        $lane_tracker++;
    }
    if ( $direction == 2 ) {
        ## 2nd read direction

        ## $lane_tracker -1 since it gets incremented after direction eq 1
        # Alias
        $mip_file_format = $infile_lane_prefix_href->{$sample_id}[ $lane_tracker - 1 ];

        my %direction_two_metric = ( sequence_run_type => q{paired-end}, );

        ## Alias
        my $file_level_href =
          \%{ $sample_info_href->{sample}{$sample_id}{file}{$mip_file_format} };
      INFO:
        while ( my ( $file_key, $file_value ) = each %direction_two_metric ) {

            $file_level_href->{$file_key} = $file_value;
        }
    }

# Save new format in hash with sample id as keys and inputfiles in array. Note: These files have not been created yet and there is one entry per strand and the file suffix is removed (.fastq).
    $infile_both_strands_prefix_href->{$sample_id}[$file_index] =
      $mip_file_format_with_direction;

    my %both_directions_metric = (
        date               => $parsed_date,
        flowcell           => $flowcell,
        lane               => $lane,
        original_file_name => $file_info_href->{$sample_id}{mip_infiles}[$file_index],
        original_file_name_prefix => $original_file_name_prefix,
        read_direction            => $direction,
        run_barcode               => $run_barcode,
        sample_barcode            => $index,
    );

    ## Alias
    my $direction_level_href =
      \%{ $sample_info_href->{sample}{$sample_id}{file}{$mip_file_format}
          {read_direction_file}{$mip_file_format_with_direction} };

  INFO:
    while ( my ( $file_key, $file_value ) = each %both_directions_metric ) {

        $direction_level_href->{$file_key} = $file_value;
    }

    return $lane_tracker;
}

sub set_most_complete_vcf {

## Function : Sets the most complete vcf file to sample_info
## Returns  :
## Arguments: $active_parameter_href     => Active parameters for this analysis hash {REF}
##          : $path                      => Path to file
##          : $recipe_name               => Recipe name
##          : $sample_info_href          => Info on samples and case hash {REF}
##          : $vcf_file_key              => Key for labelling most complete vcf
##          : $vcfparser_outfile_counter => Number of outfile files from in vcfParser (select, range)

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $path;
    my $recipe_name;
    my $sample_info_href;

    ## Default(s)
    my $vcf_file_key;
    my $vcfparser_outfile_counter;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        path => { defined => 1, required => 1, store => \$path, strict_type => 1, },
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
        vcf_file_key => {
            allow   => [qw{ vcf_file vcf_binary_file sv_vcf_file sv_vcf_binary_file }],
            default => q{vcf_file},
            store   => \$vcf_file_key,
            strict_type => 1,
        },
        vcfparser_outfile_counter => {
            default     => 0,
            store       => \$vcfparser_outfile_counter,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    if ( $active_parameter_href->{$recipe_name} == 1 ) {

        if ( $vcfparser_outfile_counter == 1 ) {

            $sample_info_href->{$vcf_file_key}{clinical}{path} = $path;
        }
        else {

            $sample_info_href->{$vcf_file_key}{research}{path} = $path;
        }
    }
    return;
}

sub set_parameter_in_sample_info {

##Function : Sets parameter to sample info
##Returns  :
##Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##         : $key_to_add            => Key and value to add
##         : $sample_info_href      => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $key_to_add;
    my $sample_info_href;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        key_to_add => {
            defined     => 1,
            required    => 1,
            store       => \$key_to_add,
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

    if ( exists $active_parameter_href->{$key_to_add} ) {

        $sample_info_href->{$key_to_add} = $active_parameter_href->{$key_to_add};
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

                $sample_info_href->{sample}{$sample_id}{recipe}{$recipe_name}
                  {$parameter_key} = $parameter_value;
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

                $sample_info_href->{recipe}{$recipe_name}{$metafile_tag}{$parameter_key}
                  = $parameter_value;
            }
        }
    }
    return;
}

sub set_in_sample_info {

## Function : Sets parameter info to sample_info
## Returns  :
## Arguments: $active_parameter_href  => Active parameters for this analysis hash {REF}
##          : $case_id_ref            => The case_id_ref {REF}
##          : $file_info_href         => File info hash {REF}
##          : $human_genome_reference => Human genome reference
##          : $outdata_dir            => Outdata directory
##          : $sample_info_href       => Info on samples and case hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $file_info_href;
    my $sample_info_href;

    ## Default(s)
    my $case_id_ref;
    my $human_genome_reference;
    my $outdata_dir;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        case_id_ref => {
            default     => \$arg_href->{active_parameter_href}{case_id},
            store       => \$case_id_ref,
            strict_type => 1,
        },
        file_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$file_info_href,
            strict_type => 1,
        },
        human_genome_reference => {
            default     => $arg_href->{active_parameter_href}{human_genome_reference},
            store       => \$human_genome_reference,
            strict_type => 1,
        },
        outdata_dir => {
            default     => $arg_href->{active_parameter_href}{outdata_dir},
            store       => \$outdata_dir,
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

    use MIP::Pedigree qw{ has_trio };

    ## Add parameter key to sample info
    my @add_keys = qw{ analysis_type expected_coverage };

  PARAMETER:
    foreach my $key_to_add (@add_keys) {

        set_parameter_in_sample_info(
            {
                active_parameter_href => $active_parameter_href,
                key_to_add            => $key_to_add,
                sample_info_href      => $sample_info_href,
            }
        );
    }

    ## Addition of genome build version to sample_info
    if ( defined $human_genome_reference ) {

        $sample_info_href->{human_genome_build}{path} = $human_genome_reference;

        my @human_genome_features = qw{ source version };
        foreach my $feature (@human_genome_features) {

            $sample_info_href->{human_genome_build}{$feature} =
              $file_info_href->{ q{human_genome_reference_} . $feature };
        }
    }
    if ( exists( $active_parameter_href->{pedigree_file} ) ) {

        ## Add pedigree_file to sample_info
        $sample_info_href->{pedigree_file}{path} =
          $active_parameter_href->{pedigree_file};
    }
    if ( exists( $active_parameter_href->{log_file} ) ) {

        ## Add log_file_dir to sample info file
        my $path = dirname( dirname( $active_parameter_href->{log_file} ) );
        $sample_info_href->{log_file_dir}       = $path;
        $sample_info_href->{last_log_file_path} = $active_parameter_href->{log_file};
    }
    ## Check for trio and set
    $sample_info_href->{has_trio} = has_trio(
        {
            active_parameter_href => $active_parameter_href,
            sample_info_href      => $sample_info_href,
        }
    );

    return;
}

sub _file_name_formats {

## Function : Define format using information derived from infile name.
## Returns  : $mip_file_format, $mip_file_format_with_direction, $original_file_name_prefix, $run_barcode;
## Arguments: $date               => Flow-cell sequencing date
##          : $direction          => Sequencing read direction
##          : $original_file_name => Original file name
##          : $flowcell           => Flow-cell id
##          : $index              => The DNA library preparation molecular barcode
##          : $lane               => Flow-cell lane
##          : $sample_id          => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $date;
    my $direction;
    my $flowcell;
    my $index;
    my $lane;
    my $original_file_name;
    my $sample_id;

    my $tmpl = {
        date => {
            defined     => 1,
            required    => 1,
            store       => \$date,
            strict_type => 1,
        },
        direction => {
            allow       => [ 1, 2 ],
            defined     => 1,
            required    => 1,
            store       => \$direction,
            strict_type => 1,
        },
        flowcell => {
            defined     => 1,
            required    => 1,
            store       => \$flowcell,
            strict_type => 1,
        },
        index => { defined => 1, required => 1, store => \$index, strict_type => 1, },
        lane  => {
            allow       => qr{ \A\d+\z }xsm,
            defined     => 1,
            required    => 1,
            store       => \$lane,
            strict_type => 1,
        },
        original_file_name => {
            defined     => 1,
            required    => 1,
            store       => \$original_file_name,
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

    my $mip_file_format =
        $sample_id
      . $UNDERSCORE
      . $date
      . $UNDERSCORE
      . $flowcell
      . $UNDERSCORE
      . $index
      . $UNDERSCORE . q{lane}
      . $lane;

    my $mip_file_format_with_direction = $mip_file_format . $UNDERSCORE . $direction;

    my $original_file_name_prefix = substr $original_file_name, 0,
      index $original_file_name, q{.fastq};

    my $run_barcode =
      $date . $UNDERSCORE . $flowcell . $UNDERSCORE . $lane . $UNDERSCORE . $index;
    return $mip_file_format, $mip_file_format_with_direction,
      $original_file_name_prefix, $run_barcode;
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
