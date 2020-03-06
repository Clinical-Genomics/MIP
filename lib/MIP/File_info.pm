package MIP::File_info;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname fileparse };
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
    our $VERSION = 1.03;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      set_alt_loci_contigs
      set_bam_contigs
      set_dict_contigs
      set_file_tag
      set_human_genome_reference_features
      set_primary_contigs
    };
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

1;
