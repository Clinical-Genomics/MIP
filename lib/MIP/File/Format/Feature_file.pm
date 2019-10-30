package MIP::File::Format::Feature_file;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
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
  qw{ $COLON $COMMA $EMPTY_STR $NEWLINE $SEMICOLON $SPACE $TAB $UNDERSCORE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ parse_feature_file_data parse_feature_file_header read_feature_file set_vcf_header_info tree_annotations };
}

sub parse_feature_file_data {

## Function : Parse feature file data and build interval tree from the data
## Returns  :
## Arguments: $data_line               => Data line
##          : $feature_columns_ref     => Feature columns to include {REF}
##          : $feature_data_href       => Feature file hash {REF}
##          : $feature_file_type       => Feature file key used to distinguish feature file(s) i.e., select or range
##          : $feature_matching_column => Column in the feature file to match with vcf key annotation {Optional}
##          : $padding                 => Padding distance
##          : $tree_href               => Interval tree hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $data_line;
    my $feature_columns_ref;
    my $feature_data_href;
    my $feature_file_type;
    my $feature_matching_column;
    my $padding;
    my $tree_href;

    my $tmpl = {
        data_line => {
            defined     => 1,
            required    => 1,
            store       => \$data_line,
            strict_type => 1,
        },
        feature_columns_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$feature_columns_ref,
            strict_type => 1,
        },
        feature_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$feature_data_href,
            strict_type => 1,
        },
        feature_file_type => {
            defined     => 1,
            required    => 1,
            store       => \$feature_file_type,
            strict_type => 1,
        },
        feature_matching_column =>
          { store => \$feature_matching_column, strict_type => 1, },
        padding => {
            defined     => 1,
            required    => 1,
            store       => \$padding,
            strict_type => 1,
        },
        tree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$tree_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Vcfparser qw{ build_interval_tree };

    ## Split data into array elements
    my @data_features = split $TAB, $data_line;

    if ( defined $feature_matching_column ) {

        my $data_feature = $data_features[$feature_matching_column];

        # Replace whitespace with underscore
        $data_feature =~ s/\s/_/gsxm;

        ## Set matching column data feature to feature data
        $feature_data_href->{$data_feature} = $data_feature;
    }

    ## Create Interval Tree
    if ( @{$feature_columns_ref} ) {

        ## Annotate vcf with features from feature file
        build_interval_tree(
            {
                feature_columns_ref => $feature_columns_ref,
                feature_file_type   => $feature_file_type,
                line_elements_ref   => \@data_features,
                padding             => $padding,
                tree_href           => $tree_href,
            }
        );
    }
    return 1;
}

sub parse_feature_file_header {

## Function : Get feature file header
## Returns  :
## Arguments: $feature_columns_ref => Feature columns to include {REF}
##          : $feature_data_href   => Feature file hash {REF}
##          : $feature_file_type   => Feature file key used to distinguish feature file(s) i.e., select or range
##          : $feature_file_path   => Feature file path
##          : $header_line         => Header line

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $feature_columns_ref;
    my $feature_data_href;
    my $feature_file_type;
    my $feature_file_path;
    my $header_line;

    my $tmpl = {
        feature_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$feature_data_href,
            strict_type => 1,
        },
        feature_columns_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$feature_columns_ref,
            strict_type => 1,
        },
        feature_file_type => {
            defined     => 1,
            required    => 1,
            store       => \$feature_file_type,
            strict_type => 1,
        },
        feature_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$feature_file_path,
            strict_type => 1,
        },
        header_line => {
            defined     => 1,
            required    => 1,
            store       => \$header_line,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Format::Feature_file qw{ set_vcf_header_info };

    ## Split headers into array elements
    my @headers = split $TAB, $header_line;

    ## Defines what headers to store from feature file
    while ( my ( $feature_index, $feature_position ) = each @{$feature_columns_ref} ) {

        ## Alias
        my $header_key = $headers[$feature_position];

        set_vcf_header_info(
            {
                feature_file_type => $feature_file_type,
                feature_file_path => $feature_file_path,
                header_key        => $header_key,
                meta_data_href    => $feature_data_href,
                position          => $feature_index,
            }
        );
    }
    return 1;
}

sub read_feature_file {

## Function : Reads a file containg features to be annotated using range queries e.g. EnsemblGeneID. Adds to metadata hash and creates an interval tree for feature.
## Returns  :
## Arguments: $feature_columns_ref     => Feature columns to include {REF}
##          : $feature_data_href       => Feature file data hash {REF}
##          : $feature_file_path       => Feature file path
##          : $feature_file_type       => Feature file type used to seperate feature file(s) e.g., select and range
##          : $feature_matching_column => Column in the feature file to match with vcf key annotation {Optional}
##          : $log                     => Log object
##          : $padding                 => Padding distance (base pairs)
##          : $tree_href               => Interval tree hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $feature_columns_ref;
    my $feature_data_href;
    my $feature_file_path;
    my $feature_file_type;
    my $feature_matching_column;
    my $log;
    my $padding;
    my $tree_href;

    my $tmpl = {
        feature_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$feature_data_href,
            strict_type => 1,
        },
        feature_columns_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$feature_columns_ref,
            strict_type => 1,
        },
        feature_file_path => {
            required    => 1,
            store       => \$feature_file_path,
            strict_type => 1,
        },
        feature_file_type => {
            defined     => 1,
            required    => 1,
            store       => \$feature_file_type,
            strict_type => 1,
        },
        feature_matching_column =>
          { store => \$feature_matching_column, strict_type => 1, },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        padding => {
            defined     => 1,
            required    => 1,
            store       => \$padding,
            strict_type => 1,
        },
        tree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$tree_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Vcfparser qw{ build_interval_tree };
    use MIP::File::Format::Feature_file
      qw{ parse_feature_file_data parse_feature_file_header };

    return if ( not $feature_file_path );

    my $filehandle = IO::Handle->new();

    open $filehandle, q{<}, $feature_file_path
      or
      $log->logdie( q{Cannot open } . $feature_file_path . $COLON . $OS_ERROR, $NEWLINE );

  LINE:
    while (<$filehandle>) {

        ## Remove newline
        chomp;

        ## Unpack line
        my $line = $_;

        ## Skip blank lines
        next LINE if ( $line =~ /^\s+$/sxm );

        ## Skip meta data lines
        next LINE if ( $line =~ /\A [#]{2}/sxm );

        ## Feature file header
        if ( $line =~ /\A [#]{1}/sxm ) {

            parse_feature_file_header(
                {
                    feature_columns_ref => $feature_columns_ref,
                    feature_data_href   => $feature_data_href,
                    feature_file_type   => $feature_file_type,
                    feature_file_path   => $feature_file_path,
                    header_line         => $line,
                }
            );
            next LINE;
        }

        ## Feature file data
        parse_feature_file_data(
            {
                data_line               => $line,
                feature_columns_ref     => $feature_columns_ref,
                feature_data_href       => $feature_data_href,
                feature_file_type       => $feature_file_type,
                padding                 => $padding,
                feature_matching_column => $feature_matching_column,
                tree_href               => $tree_href,
            }
        );
    }
    close $filehandle;
    $log->info(qq{Finished reading $feature_file_type file: $feature_file_path});
    return 1;
}

sub set_vcf_header_info {

## Function : Adds arbitrary INFO fields to hash based on supplied header key
##            unless header key is already defined
## Returns  :
## Arguments: $feature_file_type => Feature file key
##          : $feature_file_path => Feature file path
##          : $header_key        => Header key from feature file
##          : $meta_data_href    => Hash to store meta_data in {REF}
##          : $position          => Column position in supplied range file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $feature_file_type;
    my $feature_file_path;
    my $header_key;
    my $meta_data_href;
    my $position;

    my $tmpl = {
        feature_file_type => {
            defined     => 1,
            required    => 1,
            store       => \$feature_file_type,
            strict_type => 1,
        },
        feature_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$feature_file_path,
            strict_type => 1,
        },
        header_key => {
            defined     => 1,
            required    => 1,
            store       => \$header_key,
            strict_type => 1,
        },
        meta_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$meta_data_href,
            strict_type => 1,
        },
        position => {
            defined     => 1,
            required    => 1,
            store       => \$position,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## For not previously defined header keys in feature files definition
    my $arbitrary_info_field =
        q{##INFO=<ID=}
      . $header_key
      . q{,Number=.,Type=String,Description="String taken from }
      . $feature_file_path . q{">};

    ## Add INFO from predefined entries
    if ( defined $meta_data_href->{$feature_file_type}{$header_key} ) {

        $meta_data_href->{present}{$header_key}{INFO} =
          $meta_data_href->{$feature_file_type}{$header_key}{INFO};
    }
    else {
        ## Add arbitrary INFO field using feature file header key

        $meta_data_href->{present}{$header_key}{INFO} = $arbitrary_info_field;
    }

    ## Column position in supplied tsv feature file
    $meta_data_href->{present}{$header_key}{column_order} =
      $position;

    return;
}

sub tree_annotations {

## Function : Checks if an interval tree exists (per chr) and collects features
##            from input array and adds annotations to line
## Returns  :
## Arguments: $alt_allele_field  => Alternativ allele field
##          : $contig            => Contig
##          : $data_href         => Range file hash {REF}
##          : $feature_file_type => Feature file type
##          : $record_href       => Record hash info {REF}
##          : $ref_allele        => Reference allele
##          : $start             => Variant start position
##          : $tree_href         => Interval tree hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $alt_allele_field;
    my $contig;
    my $data_href;
    my $feature_file_type;
    my $record_href;
    my $ref_allele;
    my $start;
    my $tree_href;

    my $tmpl = {
        alt_allele_field => {
            defined     => 1,
            required    => 1,
            store       => \$alt_allele_field,
            strict_type => 1,
        },
        contig => {
            defined     => 1,
            required    => 1,
            store       => \$contig,
            strict_type => 1,
        },
        data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$data_href,
            strict_type => 1,
        },
        feature_file_type => {
            defined     => 1,
            required    => 1,
            store       => \$feature_file_type,
            strict_type => 1,
        },
        record_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$record_href,
            strict_type => 1,
        },
        ref_allele => {
            defined     => 1,
            required    => 1,
            store       => \$ref_allele,
            strict_type => 1,
        },
        start => {
            defined     => 1,
            required    => 1,
            store       => \$start,
            strict_type => 1,
        },
        tree_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$tree_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Format::Vcf qw{ convert_to_range };

    ## No HGNC_symbol or ensembl gene ID, but still clinically releveant e.g. non_coding_regions
    my %noid_region;

    ## Return if no feature file annotations in interval tree
    return if ( not defined $tree_href->{$feature_file_type}{$contig} );

    ## Convert precise variant to range coordinates from vcf coordinates
    my $stop = convert_to_range(
        {
            alt_allele_field => $alt_allele_field,
            reference_allele => $ref_allele,
            start_position   => $start,
        }
    );

    ## Overlapping feature elements to be collected for variant range
    my $features_ref = $tree_href->{$feature_file_type}{$contig}->fetch( $start, $stop );

    ## If no overlapping features elements returned
    return if ( not @{$features_ref} );

    ## Collect all overlapping features before adding to line
    my %collected_annotation;

  FEATURE:
    foreach my $feature ( @{$features_ref} ) {

        ## Split feature string into annotations
        my @annotations = split $SEMICOLON, $feature;

      ANNOTATION:
        while ( my ( $annotation_index, $annotation ) = each @annotations ) {

            next ANNOTATION if ( not defined $annotation );

            next ANNOTATION if ( $annotation eq $EMPTY_STR );

            push @{ $collected_annotation{$annotation_index} }, $annotation;

            ## Special case, where there is no HGNC or Ensembl gene ID
            ## but the region should be included in the select file anyway
            if ( $annotation eq q{no_hgnc_id} ) {

                ## ID = CHROM_START_REF_ALT
                my $variant_id_key = join $UNDERSCORE,
                  ( $contig, $start, $ref_allele, $alt_allele_field );
                $noid_region{$variant_id_key}++;
            }
        }
    }

  COLLECTED_ANNOTATION:
    while ( my ( $annotation_index, $annotations_ref ) = each %collected_annotation ) {

        next COLLECTED_ANNOTATION if ( not @{$annotations_ref} );

        ## All feature file annotations
      FEATURE_ANNOTATION:
        for my $feature_header ( keys %{ $data_href->{present} } ) {

            my $feature_header_col_index =
              $data_href->{present}{$feature_header}{column_order};

            ## Matching feature file header index and feature index
            next FEATURE_ANNOTATION
              if ( not $feature_header_col_index eq $annotation_index );

            ## Set annotation to vcf record
            $record_href->{ q{INFO_addition_} . $feature_file_type }{$feature_header} =
              join $COMMA, @{$annotations_ref};
        }
    }
    return %noid_region;
}

1;
