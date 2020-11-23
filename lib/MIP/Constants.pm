package MIP::Constants;

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

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.27;

    # Functions and variables which can be optionally exported

    our @EXPORT_OK = qw{
      $AMPERSAND
      %ANALYSIS
      $ASTERISK
      $AT
      $BACKTICK
      $BACKWARD_SLASH
      $CLOSE_BRACE
      $CLOSE_BRACKET
      $CLOSE_PARENTHESIS
      $COLON
      $COMMA
      @CONTAINER_BIND_PATHS
      %CONTAINER_CMD
      $CONTAINER_MANAGER
      $DASH
      $DOLLAR_SIGN
      $DOT
      $DOUBLE_QUOTE
      $EMPTY_STR
      $EQUALS
      $ESCAPE
      $FORWARD_SLASH
      $GENOME_VERSION
      $GENOME_SOURCE
      $LOG_NAME
      $MIP_VERSION
      $MOOSEX_APP_SCEEN_WIDTH
      $NEWLINE
      $OPEN_BRACE
      $OPEN_BRACKET
      $OPEN_PARENTHESIS
      $PIPE
      $PLUS
      %PRIMARY_CONTIG
      $SEMICOLON
      $SINGLE_QUOTE
      %SO_CONSEQUENCE_SEVERITY
      $SPACE
      $TAB
      $TEST_MODE
      %TEST_PROCESS_RETURN
      $UNDERSCORE
      set_container_cmd
      set_container_constants
      set_genome_build_constants
      set_log_name_constant
      set_test_constants
    };
}

## Constants

## Analysis
Readonly our %ANALYSIS => (
    ANNOTATION_DISTANCE    => 5000,
    ANNOTATION_DISTANCE_MT => 0,
    JAVA_GUEST_OS_MEMORY   => 4,
);

## Set MIP version
Readonly our $MIP_VERSION => q{v9.0.5};

## Cli
Readonly our $MOOSEX_APP_SCEEN_WIDTH => 160;

## SO-terms
Readonly our %SO_CONSEQUENCE_SEVERITY => (
    transcript_ablation => {
        rank                      => 1,
        genetic_region_annotation => q{exonic}
    },
    splice_donor_variant => {
        rank                      => 2,
        genetic_region_annotation => q{splicing}
    },
    splice_acceptor_variant => {
        rank                      => 2,
        genetic_region_annotation => q{splicing}
    },
    stop_gained => {
        rank                      => 3,
        genetic_region_annotation => q{exonic}
    },
    frameshift_variant => {
        rank                      => 4,
        genetic_region_annotation => q{exonic}
    },
    stop_lost => {
        rank                      => 5,
        genetic_region_annotation => q{exonic}
    },
    start_lost => {
        rank                      => 5,
        genetic_region_annotation => q{exonic}
    },
    initiator_codon_variant => {
        rank                      => 6,
        genetic_region_annotation => q{exonic}
    },
    inframe_insertion => {
        rank                      => 6,
        genetic_region_annotation => q{exonic}
    },
    inframe_deletion => {
        rank                      => 6,
        genetic_region_annotation => q{exonic}
    },
    missense_variant => {
        rank                      => 6,
        genetic_region_annotation => q{exonic}
    },
    protein_altering_variant => {
        rank                      => 6,
        genetic_region_annotation => q{exonic}
    },
    transcript_amplification => {
        rank                      => 7,
        genetic_region_annotation => q{exonic}
    },
    splice_region_variant => {
        rank                      => 8,
        genetic_region_annotation => q{splicing}
    },
    incomplete_terminal_codon_variant => {
        rank                      => 9,
        genetic_region_annotation => q{exonic}
    },
    synonymous_variant => {
        rank                      => 10,
        genetic_region_annotation => q{exonic}
    },
    stop_retained_variant => {
        rank                      => 10,
        genetic_region_annotation => q{exonic}
    },
    start_retained_variant => {
        rank                      => 10,
        genetic_region_annotation => q{exonic}
    },
    coding_sequence_variant => {
        rank                      => 11,
        genetic_region_annotation => q{exonic}
    },
    mature_miRNA_variant => {
        rank                      => 12,
        genetic_region_annotation => q{ncRNA_exonic}
    },
    q{5_prime_UTR_variant} => {
        rank                      => 13,
        genetic_region_annotation => q{5UTR}
    },
    q{3_prime_UTR_variant} => {
        rank                      => 14,
        genetic_region_annotation => q{3UTR}
    },
    non_coding_transcript_exon_variant => {
        rank                      => 15,
        genetic_region_annotation => q{ncRNA_exonic}
    },
    non_coding_transcript_variant => {
        rank                      => 15,
        genetic_region_annotation => q{ncRNA}
    },
    intron_variant => {
        rank                      => 16,
        genetic_region_annotation => q{intronic}
    },
    NMD_transcript_variant => {
        rank                      => 17,
        genetic_region_annotation => q{ncRNA}
    },
    upstream_gene_variant => {
        rank                      => 18,
        genetic_region_annotation => q{upstream}
    },
    downstream_gene_variant => {
        rank                      => 19,
        genetic_region_annotation => q{downstream}
    },
    TFBS_ablation => {
        rank                      => 20,
        genetic_region_annotation => q{TFBS}
    },
    TFBS_amplification => {
        rank                      => 21,
        genetic_region_annotation => q{TFBS}
    },
    TF_binding_site_variant => {
        rank                      => 22,
        genetic_region_annotation => q{TFBS}
    },
    regulatory_region_variant => {
        rank                      => 22,
        genetic_region_annotation => q{regulatory_region}
    },
    regulatory_region_ablation => {
        rank                      => 23,
        genetic_region_annotation => q{regulatory_region}
    },
    regulatory_region_amplification => {
        rank                      => 24,
        genetic_region_annotation => q{regulatory_region}
    },
    feature_elongation => {
        rank                      => 25,
        genetic_region_annotation => q{genomic_feature}
    },
    feature_truncation => {
        rank                      => 26,
        genetic_region_annotation => q{genomic_feature}
    },
    intergenic_variant => {
        rank                      => 27,
        genetic_region_annotation => q{intergenic}
    },
);

## Contigs
Readonly our %PRIMARY_CONTIG => (
    38 => {
        contigs => [
            qw{
              chr1 chr2 chr3 chr4 chr5 chr6
              chr7 chr8 chr9 chr10 chr11 chr12
              chr13 chr14 chr15 chr16 chr17 chr18
              chr19 chr20 chr21 chr22 chrX chrY
              chrM },
        ],
        contigs_size_ordered => [
            qw{
              chr1 chr2 chr3 chr4 chr5 chr6
              chr7 chrX chr8 chr9 chr10 chr11
              chr12 chr13 chr14 chr15 chr16 chr17
              chr18 chr19 chr20 chr21 chr22 chrY
              chrM }
        ],
        synonyms_map => {
            chr1  => 1,
            chr2  => 2,
            chr3  => 3,
            chr4  => 4,
            chr5  => 5,
            chr6  => 6,
            chr7  => 7,
            chr8  => 8,
            chr9  => 9,
            chr10 => 10,
            chr11 => 11,
            chr12 => 12,
            chr13 => 13,
            chr14 => 14,
            chr15 => 15,
            chr16 => 16,
            chr17 => 17,
            chr18 => 18,
            chr19 => 19,
            chr20 => 20,
            chr21 => 21,
            chr22 => 22,
            chrX  => q{X},
            chrY  => q{Y},
            chrM  => q{MT},
        },
    },
    37 => {
        contigs => [
            qw{
              1 2 3 4 5 6 7 8 9 10
              11 12 13 14 15 16 17 18 19 20
              21 22 X Y MT },
        ],
        contigs_size_ordered => [
            qw{
              1 2 3 4 5 6 7 X 8 9
              10 11 12 13 14 15 16 17 18 19
              20 21 22 Y MT },
        ],
    },
    19 => {
        contigs => [
            qw{
              chr1 chr2 chr3 chr4 chr5 chr6
              chr7 chr8 chr9 chr10 chr11 chr12
              chr13 chr14 chr15 chr16 chr17 chr18
              chr19 chr20 chr21 chr22 chrX chrY
              chrM },
        ],
        contigs_size_ordered => [
            qw{
              chr1 chr2 chr3 chr4 chr5 chr6
              chr7 chrX chr8 chr9 chr10 chr11
              chr12 chr13 chr14 chr15 chr16 chr17
              chr18 chr19 chr20 chr21 chr22 chrY
              chrM }
        ],
    },
);

## Symbols
Readonly our $AMPERSAND         => q{&};
Readonly our $ASTERISK          => q{*};
Readonly our $AT                => q{@};
Readonly our $BACKTICK          => q{`};
Readonly our $BACKWARD_SLASH    => q{\\};
Readonly our $CLOSE_BRACE       => q{\}};
Readonly our $CLOSE_BRACKET     => q{]};
Readonly our $CLOSE_PARENTHESIS => q{)};
Readonly our $COLON             => q{:};
Readonly our $COMMA             => q{,};
Readonly our $DASH              => q{-};
Readonly our $DOLLAR_SIGN       => q{$};
Readonly our $DOT               => q{.};
Readonly our $DOUBLE_QUOTE      => q{"};
Readonly our $EMPTY_STR         => q{};
Readonly our $ESCAPE            => q{\\};
Readonly our $EQUALS            => q{=};
Readonly our $FORWARD_SLASH     => q{/};
Readonly our $NEWLINE           => qq{\n};
Readonly our $OPEN_BRACE        => q{\{};
Readonly our $OPEN_BRACKET      => q{[};
Readonly our $OPEN_PARENTHESIS  => q{(};
Readonly our $PIPE              => q{|};
Readonly our $PLUS              => q{+};
Readonly our $SEMICOLON         => q{;};
Readonly our $SINGLE_QUOTE      => q{'};
Readonly our $SPACE             => q{ };
Readonly our $TAB               => qq{\t};
Readonly our $UNDERSCORE        => q{_};

sub set_container_cmd {

## Function : Set container commands for each executable
## Returns  :
## Arguments: $container_cmd_href => Analysis recipe hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $container_cmd_href;

    my $tmpl = {
        container_cmd_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$container_cmd_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    Readonly our %CONTAINER_CMD => %{$container_cmd_href};

    return;
}

sub set_container_constants {

## Function : Set container constants
## Returns  :
## Arguments: $active_parameter_href => Analysis recipe hash {REF}

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

    use MIP::Environment::Path qw{ reduce_dir_paths };
    use MIP::Test::Fixtures qw{ test_constants };

    if ( $active_parameter_href->{test_mode} ) {

        test_constants( {} );
    }

    my @container_bind_paths = (
        $active_parameter_href->{reference_dir},
        $active_parameter_href->{outdata_dir},
        $active_parameter_href->{temp_directory},
    );

    @container_bind_paths = reduce_dir_paths(
        {
            dir_paths_ref => \@container_bind_paths,
        }
    );

    Readonly our @CONTAINER_BIND_PATHS => @container_bind_paths;
    Readonly our $CONTAINER_MANAGER    => $active_parameter_href->{container_manager};

    return;
}

sub set_genome_build_constants {

## Function : Set human build constants for analysis run
## Returns  :
## Arguments: $genome_source  => Source of genome (grch or hg)
##          : $genome_version => Genome build version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $genome_source;
    my $genome_version;

    my $tmpl = {
        genome_source => {
            defined     => 1,
            required    => 1,
            store       => \$genome_source,
            strict_type => 1,
        },
        genome_version => {
            defined     => 1,
            required    => 1,
            store       => \$genome_version,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Set global genome version and source
    Readonly our $GENOME_VERSION => $genome_version;
    Readonly our $GENOME_SOURCE  => $genome_source;

    return;
}

sub set_log_name_constant {

## Function : Set log name constant for MIP
## Returns  :
## Arguments: $log_name => Name of log

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $log_name;

    my $tmpl = {
        log_name => {
            defined     => 1,
            required    => 1,
            store       => \$log_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Log name
    Readonly our $LOG_NAME => $log_name;

    return;
}

sub set_test_constants {

    ## Function : Set constants for testing
    ## Returns  :
    ## Arguments: $container_manager        => Name of container manager
    ##          : $test_mode                => Set test mode
    ##          : $test_process_return_href => Process return hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $test_process_return_href;

    ## Default(s)
    my $container_manager;
    my $test_mode;

    my $tmpl = {
        container_manager => {
            allow       => [qw{docker singularity}],
            default     => q{docker},
            store       => \$container_manager,
            strict_type => 1,
        },
        test_mode => {
            allow       => [ 0, 1 ],
            default     => 1,
            store       => \$test_mode,
            strict_type => 1,
        },
        test_process_return_href => {
            default     => {},
            store       => \$test_process_return_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    our $CONTAINER_MANAGER   = $container_manager;
    our $TEST_MODE           = $test_mode;
    our %TEST_PROCESS_RETURN = %{$test_process_return_href};

    return;
}
1;
