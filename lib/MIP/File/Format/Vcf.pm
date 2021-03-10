package MIP::File::Format::Vcf;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $COLON $COMMA $DOT $EQUALS $PIPE $SEMICOLON $SPACE $TAB };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      check_vcf_variant_line
      convert_to_range
      get_transcript_effects
      parse_vcf_header
      get_bcftools_norm_command_from_vcf_header
      get_vcf_header_line_by_id
      set_in_consequence_hash
      set_info_key_pairs_in_vcf_record
      set_line_elements_in_vcf_record
    };
}

## Constants
Readonly my $INFO_COL_NR => 7;

sub get_bcftools_norm_command_from_vcf_header {

## Function : Get vcf header line matching bcftools norm command
## Returns  :
## Arguments: $vcf_file_path        => Path to vcf file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $vcf_file_path;

    my $tmpl = {
        vcf_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$vcf_file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Environment::Child_process qw{ child_process };
    use MIP::Environment::Executable qw{ get_executable_base_command };
    use MIP::Language::Perl qw{ perl_nae_oneliners };

    my @check_header_cmds;
    my $bcftools_binary_path = get_executable_base_command( { base_command => q{bcftools}, } );

    ## Stream vcf using bcftools
    push @check_header_cmds, $bcftools_binary_path . $SPACE . q{view} . $SPACE . $vcf_file_path;
    push @check_header_cmds, $PIPE;

    ## Assemble perl regexp for detecting keys in vcf
    push @check_header_cmds,
      perl_nae_oneliners(
        {
            oneliner_name      => q{bcftools_norm_check},
            oneliner_parameter => q{bcftools_norm},
        }
      );
    push @check_header_cmds, $SEMICOLON;

    my %process_return = child_process(
        {
            commands_ref => \@check_header_cmds,
            process_type => q{open3},
        }
    );

    return $process_return{stdouts_ref}[0];
}

sub get_vcf_header_line_by_id {

## Function : Get vcf header line matching given ID
## Returns  :
## Arguments: $bcftools_binary_path => Path to bcftools
##          : $header_id            => Header ID
##          : $vcf_file_path        => Path to vcf file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $bcftools_binary_path;
    my $header_id;
    my $vcf_file_path;

    my $tmpl = {
        bcftools_binary_path => {
            defined     => 1,
            required    => 1,
            store       => \$bcftools_binary_path,
            strict_type => 1,
        },
        header_id => {
            defined     => 1,
            required    => 1,
            store       => \$header_id,
            strict_type => 1,
        },
        vcf_file_path => {
            defined     => 1,
            required    => 1,
            store       => \$vcf_file_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Environment::Child_process qw{ child_process };
    use MIP::Language::Perl qw{ perl_nae_oneliners };

    my @check_header_cmds;

    ## Stream vcf using bcftools
    push @check_header_cmds, $bcftools_binary_path . $SPACE . q{view} . $SPACE . $vcf_file_path;
    push @check_header_cmds, $PIPE;

    ## Assemble perl regexp for detecting keys in vcf
    push @check_header_cmds,
      perl_nae_oneliners(
        {
            oneliner_name      => q{get_vcf_header_id_line},
            oneliner_parameter => $header_id,
        }
      );
    push @check_header_cmds, $SEMICOLON;

    my %process_return = child_process(
        {
            commands_ref => \@check_header_cmds,
            process_type => q{open3},
        }
    );

    return $process_return{stdouts_ref}[0];
}

sub check_vcf_variant_line {

## Function : Check variant line elements
## Returns  :
## Arguments: $input_line_number         => Input line number
##          : $log                       => Log object
##          : $variant_line              => Variant line
##          : $variant_line_elements_ref => Array for variant line elements {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $input_line_number;
    my $log;
    my $variant_line;
    my $variant_line_elements_ref;

    my $tmpl = {
        input_line_number => {
            defined     => 1,
            required    => 1,
            store       => \$input_line_number,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        variant_line => {
            defined     => 1,
            required    => 1,
            store       => \$variant_line,
            strict_type => 1,
        },
        variant_line_elements_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$variant_line_elements_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Check that we have some data in INFO field
    return 1 if ( defined $variant_line_elements_ref->[$INFO_COL_NR] );

    $log->fatal(qq{No INFO field at line number: $input_line_number});
    $log->fatal(qq{Displaying malformed line: $variant_line});
    exit 1;
}

sub convert_to_range {

## Function : Converts VCF variants to corresponding range coordinates
## Returns  : $final_stop_position
## Arguments: $alt_allele_field => Alternative allele field
##          : $reference_allele => Reference allele
##          : $start_position   => Variant start position

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $alt_allele_field;
    my $reference_allele;
    my $start_position;

    my $tmpl = {
        alt_allele_field => {
            defined     => 1,
            required    => 1,
            store       => \$alt_allele_field,
            strict_type => 1,
        },
        reference_allele => {
            defined     => 1,
            required    => 1,
            store       => \$reference_allele,
            strict_type => 1,
        },
        start_position => {
            defined     => 1,
            required    => 1,
            store       => \$start_position,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## No alternative allele call
    return if ( $alt_allele_field eq $DOT );

    ## The most "downstream" position per variant
    my $final_stop_position = 0;

    ## Split into alternative allele(s)
    my @alt_alleles = split $COMMA, $alt_allele_field;

  ALLELE:
    foreach my $alt_allele (@alt_alleles) {

        my $stop_position;

        ## SNV
        if (    length $reference_allele == 1
            and length $alt_allele == 1 )
        {

            $stop_position = $start_position + 1;
        }
        elsif ( length $reference_allele >= length $alt_allele ) {
            ## Deletion or block substitution
            $stop_position = $start_position + length($reference_allele) - 1;
        }
        elsif ( length $reference_allele < length $alt_allele ) {
            ## insertion or block substitution
            $stop_position = $start_position + length($alt_allele) - 1;
        }

        ## Collect largest range per variant record based on all alternative_alleles
        # New end is downstream of old
        if ( $final_stop_position < $stop_position ) {

            $final_stop_position = $stop_position;
        }
    }
    return $final_stop_position;
}

sub get_transcript_effects {

## Function : Set and return csq hash per transcript
## Returns  :
## Arguments: $transcript_effects_ref       => Transcript effects {REF}
##          : $vep_format_field_column_href => VEP format columns {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $transcript_effects_ref;
    my $vep_format_field_column_href;

    my $tmpl = {
        transcript_effects_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$transcript_effects_ref,
            strict_type => 1,
        },
        vep_format_field_column_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$vep_format_field_column_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Vcfparser qw{ %CSQ_FIELD_MAP };

    my %csq;

  CSQ_FIELD:
    while ( my ( $csq_format_field_key, $csq_key_id ) = each %CSQ_FIELD_MAP ) {

        ## Get where CSQ annotation is in transcript_effects according to VEP CSQ Schema
        my $annotation_index = $vep_format_field_column_href->{$csq_format_field_key};
        $csq{$csq_key_id} = $transcript_effects_ref->[$annotation_index];
    }
    return %csq;
}

sub parse_vcf_header {

## Function : Adds header meta data to hash
## Returns  :
## Arguments: $meta_data_href   => Hash for header data {REF}
##          : $meta_data_string => Meta data string from vcf header

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $meta_data_href;
    my $meta_data_string;

    my $tmpl = {
        meta_data_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$meta_data_href,
            strict_type => 1,
        },
        meta_data_string => {
            defined     => 1,
            required    => 1,
            store       => \$meta_data_string,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Define how to parse meta-data header
    my %vcf_header_regexp = (
        fileformat   => q{\A [#]{2}(fileformat)=(\S+)},
        field_format => q{\A [#]{2}([^=]*)=<ID=([^,]*)},
    );

  SCHEMA:
    while ( my ( $key, $regexp ) = each %vcf_header_regexp ) {

        my ( $vcf_schema, $vcf_id ) = $meta_data_string =~ /$regexp/sxm;

        ## Alias method subs
        my $set_by_schema_and_id_cref = sub {
            $meta_data_href->{$vcf_schema}{$vcf_id} = $meta_data_string;
        };
        my $add_by_schema_cref = sub {
            push @{ $meta_data_href->{$vcf_schema}{$vcf_schema} }, $meta_data_string;
        };

        ## Create dispatch table
        my %add_to_meta_data = (
            ALT        => $set_by_schema_and_id_cref,
            contig     => $add_by_schema_cref,
            fileformat => $set_by_schema_and_id_cref,
            FILTER     => $set_by_schema_and_id_cref,
            FORMAT     => $set_by_schema_and_id_cref,
            INFO       => $set_by_schema_and_id_cref,
        );

        if ( $vcf_schema and exists $add_to_meta_data{$vcf_schema} ) {

            ## Add header line to hash
            $add_to_meta_data{$vcf_schema}->();
            return;
        }
    }

    ## All other meta-data headers - Add to array without using regexp
    push @{ $meta_data_href->{other}{other} }, $meta_data_string;
    return 1;
}

sub set_in_consequence_hash {

## Function : Set most severe consequence key set in hash
## Returns  :
## Arguments: $allele           => Allele
##          : $consequence_href => Consequence hash {REF}
##          : $hgnc_id          => Hgnc id
##          : $set_key_href     => Key value pairs to set {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $allele;
    my $consequence_href;
    my $hgnc_id;
    my $set_key_href;

    my $tmpl = {
        allele => {
            defined     => 1,
            required    => 1,
            store       => \$allele,
            strict_type => 1,
        },
        consequence_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$consequence_href,
            strict_type => 1,
        },
        hgnc_id => {
            defined     => 1,
            required    => 1,
            store       => \$hgnc_id,
            strict_type => 1,
        },
        set_key_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$set_key_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

  KEY_VALUE_PAIR:
    while ( my ( $key, $value ) = each %{$set_key_href} ) {

        $consequence_href->{$hgnc_id}{$allele}{$key} = $value;
    }
    return;
}

sub set_info_key_pairs_in_vcf_record {

## Function : Adds INFO key value pairs to record hash
## Returns  :
## Arguments: $vcf_record_href => Hash for variant line data {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $vcf_record_href;

    my $tmpl = {
        vcf_record_href => {
            default  => {},
            defined  => 1,
            required => 1,
            store    => \$vcf_record_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Add INFO elements
    my @info_elements = split $SEMICOLON, $vcf_record_href->{INFO};

    ## Collect key value pairs from INFO field elements
  ELEMENT:
    foreach my $element (@info_elements) {

        my ( $key, $value ) = split $EQUALS, $element;

        $vcf_record_href->{INFO_key_value}{$key} = $value;
    }
    return;
}

sub set_line_elements_in_vcf_record {

## Function : Adds variant line elements to record hash
## Returns  :
## Arguments: $line_elements_ref      => Variant line elements {REF}
##          : $vcf_format_columns_ref => VCF format colums
##          : $vcf_record_href        => Hash for variant line data {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $line_elements_ref;
    my $vcf_format_columns_ref;
    my $vcf_record_href;

    my $tmpl = {
        line_elements_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$line_elements_ref,
            strict_type => 1,
        },
        vcf_format_columns_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$vcf_format_columns_ref,
            strict_type => 1,
        },
        vcf_record_href => {
            default  => {},
            defined  => 1,
            required => 1,
            store    => \$vcf_record_href,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Add vcf_schema as keys and line elements as value to record hash
    @{$vcf_record_href}{ @{$vcf_format_columns_ref} } = @{$line_elements_ref};

    return;
}

1;
