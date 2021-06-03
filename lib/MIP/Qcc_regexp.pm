package MIP::Qcc_regexp;

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
use MIP::Constants qw{ $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ get_qcc_regexp_recipe_attribute regexp_to_yaml };
}

sub get_qcc_regexp_recipe_attribute {

## Function : Get recipe attributes from qccollect reg exp hash
## Returns  : "$attribute" or "$attribute_href"
## Arguments: $attribute       => Attribute key
##          : $recipe_name     => Recipe to get attributes from
##          : $qcc_regexp_href => qccollect regexp hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $attribute;
    my $recipe_name;
    my $qcc_regexp_href;

    my $tmpl = {
        attribute => {
            store       => \$attribute,
            strict_type => 1,
        },
        recipe_name => {
            defined     => 1,
            required    => 1,
            store       => \$recipe_name,
            strict_type => 1,
        },
        qcc_regexp_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$qcc_regexp_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Get and return attribute value
    if ( defined $attribute && $attribute ) {

        return $qcc_regexp_href->{$recipe_name}{$attribute};
    }

    ## Get recipe attribute hash
    return %{ $qcc_regexp_href->{$recipe_name} };
}

sub regexp_to_yaml {

## Function : Write default regexp to YAML
## Returns  :
## Arguments: $log                  => Log
##          : $print_regexp_outfile => File to print regexp to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $log;
    my $print_regexp_outfile;

    my $tmpl = {
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        print_regexp_outfile => {
            required    => 1,
            store       => \$print_regexp_outfile,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Io::Write qw{ write_to_file };

    return 1 if ( not defined $print_regexp_outfile );

    my %regexp;

    ## Add to %regexp to enable print in YAML

    # Return Uniquely mapped reads %
    $regexp{star_log}{percentage_uniquely_mapped_reads} =
      q?perl -nae 'if(m/Uniquely\smapped\sreads\s%\s\|\t(\d+\.\d+) /xms) {print $1; last;}' ?;

    ## Return percentage of reads with adapters
    $regexp{trim_galore_stats}{percentage_reads_with_adapter} =
      q?perl -nae 'if( m/Reads\swith\sadapters[^(]+\((\d+\.\d+) /xms ){ print $1; last;}' ?;

    ## Return percentage of reads after trimming
    $regexp{trim_galore_stats}{percentage_reads_after_trimming} =
      q?perl -nae 'if( m/Reads\swritten\s\([^(]+\((\d+\.\d+) /xms ){ print $1; last;}' ?;

    ## Return percentage of bp remaining after trimming
    $regexp{trim_galore_stats}{percentage_bp_after_trimming} =
      q?perl -nae 'if( m/Total\swritten\s\([^(]+\((\d+\.\d+) /xms ){ print $1; last;}' ?;

    # Return Encoding
    $regexp{fastqc_ar}{encoding} =
q?perl -nae' if ($_=~/Encoding\s+(\S+\s\S+\s\S+\s\S+|\S+\s\S+)/) { my $encoding = $1;$encoding=~s/\s/\_/g; print $encoding;last;}' ?;

    # Return Sequence length
    $regexp{fastqc_ar}{sequence_length} =
      q?perl -nae' if ($_=~/Sequence length\s(\d+)/) {print $1;last;}' ?;

    # Return Total sequences
    $regexp{fastqc_ar}{total_number_of_reads} =
      q?perl -nae' if ($_=~/Total Sequences\s(\d+)/) {print $1;last;}' ?;

    # Return GC content
    $regexp{fastqc_ar}{gc} = q?perl -nae' if ($_=~/%GC\s(\d+)/) {print $1;last;}' ?;

    # Return Sequence duplication level
    $regexp{fastqc_ar}{sequence_duplication} =
      q?perl -nae' if ($_=~/#Total Duplicate Percentage\s+(\d+.\d)/) {print $1;last;}' ?;

    # Return Basic Statistics
    $regexp{fastqc_ar}{basic_statistics} =
      q?perl -nae' if ($_=~/>>Basic Statistics\s+(\S+)/) {print $1;last;}' ?;

    # Return Per base sequence quality
    $regexp{fastqc_ar}{per_base_sequence_quality} =
      q?perl -nae' if ($_=~/>>Per base sequence quality\s+(\S+)/) {print $1;last;}' ?;

    # Return Per sequence quality scores
    $regexp{fastqc_ar}{per_sequence_quality_scores} =
      q?perl -nae' if ($_=~/>>Per sequence quality scores\s+(\S+)/) {print $1;last;}' ?;

    # Return Per base sequence content
    $regexp{fastqc_ar}{per_base_sequence_content} =
      q?perl -nae' if ($_=~/>>Per base sequence content\s+(\S+)/) {print $1;last;}' ?;

    # Return Per base GC content
    $regexp{fastqc_ar}{per_base_gc_content} =
      q?perl -nae' if ($_=~/>>Per base GC content\s+(\S+)/) {print $1;last;}' ?;

    # Return Per sequence GC content
    $regexp{fastqc_ar}{per_sequence_gc_content} =
      q?perl -nae' if ($_=~/>>Per sequence GC content\s+(\S+)/) {print $1;last;}' ?;

    # Return Per base N content
    $regexp{fastqc_ar}{per_base_n_content} =
      q?perl -nae' if ($_=~/>>Per base N content\s+(\S+)/) {print $1;last;}' ?;

    # Return Sequence Duplication Levels
    $regexp{fastqc_ar}{sequence_duplication_levels} =
      q?perl -nae' if ($_=~/>>Sequence Duplication Levels\s+(\S+)/) {print $1;last;}' ?;

    # Return Overrepresented sequences
    $regexp{fastqc_ar}{overrepresented_sequences} =
      q?perl -nae' if ($_=~/>>Overrepresented sequences\s+(\S+)/) {print $1;last;}' ?;

    # Return Kmer Content
    $regexp{fastqc_ar}{kmer_content} =
      q?perl -nae' if ($_=~/>>Kmer Content\s+(\S+)/) {print $1;last;}' ?;

    # Return % mapped reads from BAM alignment
    $regexp{bamstats}{percentage_mapped_reads} =
      q?perl -nae 'if($_=~/percentage mapped reads:\s+(\S+)/) {print $1;last}' ?;

    # Return raw total sequences from BAM alignment
    $regexp{bamstats}{raw_total_sequences} =
      q?perl -nae 'if($_=~/raw total sequences:\s+(\S+)/) {print $1;last}' ?;

    # Return reads mapped from BAM alignment
    $regexp{bamstats}{reads_mapped} =
      q?perl -nae 'if($_=~/reads mapped:\s+(\S+)/) {print $1;last}' ?;

    # Return gender from chanjo_sexcheck
    $regexp{chanjo_sexcheck}{gender} =
      q?perl -nae 'if( ($F[0]!~/^#/) && ($F[2] =~/\S+/) ) {print $F[2];}' ?;

    # Return sample order from vcf file used to create ".ped", ".map" and hence ".mibs".
    $regexp{pedigree_check}{sample_order} =
q?perl -nae 'if ($_=~/^#CHROM/) {chomp $_; my @line = split(/\t/,$_); for (my $sample=9;$sample<scalar(@line);$sample++) { print $line[$sample], "\t";}last;}' ?;

    $regexp{inbreeding_factor}{sample_inbreeding_factor} =
q?perl -nae 'my @inbreedingFactor; if ($. > 1) {my @temp = split(/\s/,$_);push(@inbreedingFactor, $F[0].":".$F[5]); print $inbreedingFactor[0], "\t"; }' ?;

    $regexp{plink_sexcheck}{sample_sexcheck} =
q?perl -nae 'my @sexCheckFactor; if ($. > 1) {my @temp = split(/\s+/,$_);push(@sexCheckFactor,$temp[2].":".$temp[4]); print $sexCheckFactor[0], "\t"; }' ?;

    # Get entire sample relation check file
    $regexp{relation_check}{sample_relation_check} = q?perl -nae 'print $_;' ?;

    # Return fraction duplicates
    $regexp{markduplicates}{fraction_duplicates} =
      q?perl -nae 'if($_=~/Fraction Duplicates\: (\S+)/) {print $1;}' ?;

    # Get BAIT_SET line from header
    $regexp{collecthsmetrics}{header} = q?perl -nae' if ($_ =~/^BAIT_SET/ ) {print $_;last;}' ?;

    # Return line and only look at line 8 in file, where the data action is
    $regexp{collecthsmetrics}{data} =
      q?perl -nae' if ( ($. ==8) && ($_ =~/(\S+)/) ) {print $_;last;}' ?;

    # Return CATEGORY line from header
    $regexp{collectmultiplemetrics}{header} =
      q?perl -nae' if ($_ =~/^CATEGORY/ ) {print $_;last;}' ?;

    # Return FIRST_OF_PAIR
    $regexp{collectmultiplemetrics}{first_of_pair} =
      q?perl -nae' if ($_ =~/^FIRST_OF_PAIR/ ) {print $_;last;}' ?;

    # Return SECOND_OF_PAIR
    $regexp{collectmultiplemetrics}{second_of_pair} =
      q?perl -nae' if ($_ =~/^SECOND_OF_PAIR/ ) {print $_;last;}' ?;

    # Return PAIR line
    $regexp{collectmultiplemetrics}{pair} = q?perl -nae' if ($_ =~/^PAIR/ ) {print $_;last;}'  ?;

    # Return MEDIAN_INSERT_SIZE line from header
    $regexp{collectmultiplemetricsinsertsize}{header} =
      q?perl -nae' if ($_ =~/^MEDIAN_INSERT_SIZE/ ) {print $_;last;}' ?;

    # Return line and only look at line 8 in file, where the data action is
    $regexp{collectmultiplemetricsinsertsize}{data} =
      q?perl -nae' if ( ($. ==8) && ($_ =~/(\S+)/) ) {print $_;last;}' ?;

    # Get PF_BASES line from header
    $regexp{collectrnaseqmetrics}{header} = q?perl -nae' if ($_ =~/^PF_BASES/ ) {print $_;last;}' ?;

    # Return line and only look at line 8 in file, where the data action is
    $regexp{collectrnaseqmetrics}{data} =
      q?perl -nae' if ( ($. ==8) && ($_ =~/(\S+)/) ) {print $_;last;}' ?;

# Return CompOverlap CompFeatureInput line and only look at line 8, where the data action is in header
    $regexp{variantevalall}{comp_overlap_data_header} =
      q?perl -nae' if ($_ =~/^CompOverlap\s+CompFeatureInput/ ) {print $_;last;}' ?;

    # Return CompOverlap and all and none line
    $regexp{variantevalall}{comp_overlap_data_all} =
      q?perl -nae' if ( ($_ =~/^CompOverlap/) && ($_ =~/all/) && ($_ =~/none/)) {print $_;last;}' ?;

    # Return CompOverlap and known line
    $regexp{variantevalall}{comp_overlap_data_known} =
      q?perl -nae' if ( ($_ =~/^CompOverlap/) && ($_ =~/known\s/) ) {print $_;last;}' ?;

    # Return CompOverlap and novel line
    $regexp{variantevalall}{comp_overlap_data_novel} =
      q?perl -nae' if ( ($_ =~/^CompOverlap/) && ($_ =~/novel\s/) ) {print $_;last;}' ?;

    # Return CountVariants and CompFeatureInput line from header
    $regexp{variantevalall}{count_variants_data_header} =
      q?perl -nae' if ($_ =~/^CountVariants\s+CompFeatureInput/ ) {print $_;last;}' ?;

    # Return CountVariants and all line
    $regexp{variantevalall}{count_variants_data_all} =
      q?perl -nae' if ( ($_ =~/^CountVariants/) && ($_ =~/all\s/) ) {print $_;last;}' ?;

    # Return CountVariants and known line
    $regexp{variantevalall}{count_variants_data_known} =
      q?perl -nae' if ( ($_ =~/^CountVariants/) && ($_ =~/known\s/) ) {print $_;last;}' ?;

    # Return CountVariants and novel line
    $regexp{variantevalall}{count_variants_data_novel} =
      q?perl -nae' if ( ($_ =~/^CountVariants/) && ($_ =~/novel\s/) ) {print $_;last;}' ?;

    # Return IndelSummary and CompFeatureInput line from header
    $regexp{variantevalall}{indel_summary_data_header} =
      q?perl -nae' if ($_ =~/^IndelSummary\s+CompFeatureInput/ ) {print $_;last;}' ?;

    # Return IndelSummary and all line
    $regexp{variantevalall}{indel_summary_data_all} =
      q?perl -nae' if ( ($_ =~/^IndelSummary/) && ($_ =~/all\s/) ) {print $_;last;}' ?;

    # Return IndelSummary and known line
    $regexp{variantevalall}{indel_summary_data_known} =
      q?perl -nae' if ( ($_ =~/^IndelSummary/) && ($_ =~/known\s/) ) {print $_;last;}' ?;

    # Return IndelSummary and novel line
    $regexp{variantevalall}{indel_summary_data_novel} =
      q?perl -nae' if ( ($_ =~/^IndelSummary/) && ($_ =~/novel\s/) ) {print $_;last;}' ?;

    # Return MultiallelicSummary and CompFeatureInput line from header
    $regexp{variantevalall}{multiallelic_summary_data_header} =
      q?perl -nae' if ($_ =~/^MultiallelicSummary\s+CompFeatureInput/ ) {print $_;last;}' ?;

    # Return MultiallelicSummary and all line
    $regexp{variantevalall}{multiallelic_summary_data_all} =
      q?perl -nae' if ( ($_ =~/^MultiallelicSummary/) && ($_ =~/all\s/) ) {print $_;last;}' ?;

    # Return MultiallelicSummary and known line
    $regexp{variantevalall}{multiallelic_summary_data_known} =
      q?perl -nae' if ( ($_ =~/^MultiallelicSummary/) && ($_ =~/known\s/) ) {print $_;last;}' ?;

    # Return MultiallelicSummary and novel line
    $regexp{variantevalall}{multiallelic_summary_data_novel} =
      q?perl -nae' if ( ($_ =~/^MultiallelicSummary/) && ($_ =~/novel\s/) ) {print $_;last;}' ?;

    # Return TiTvVariantEvaluator and CompFeatureInput line from header
    $regexp{variantevalall}{titv_variant_evaluator_data_header} =
      q?perl -nae' if ($_ =~/^TiTvVariantEvaluator\s+CompFeatureInput/ ) {print $_;last;}' ?;

    # Return TiTvVariantEvaluator and all line
    $regexp{variantevalall}{titv_variant_evaluator_data_all} =
      q?perl -nae' if ( ($_ =~/^TiTvVariantEvaluator/) && ($_ =~/all\s/) ) {print $_;last;}' ?;

    # Return TiTvVariantEvaluator and known line
    $regexp{variantevalall}{titv_variant_evaluator_data_known} =
      q?perl -nae' if ( ($_ =~/^TiTvVariantEvaluator/) && ($_ =~/known\s/) ) {print $_;last;}' ?;

    # Return TiTvVariantEvaluator and novel line
    $regexp{variantevalall}{titv_variant_evaluator_data_novel} =
      q?perl -nae' if ( ($_ =~/^TiTvVariantEvaluator/) && ($_ =~/novel\s/) ) {print $_;last;}' ?;

    # Return ValidationReport and CompFeatureInput line from header
    $regexp{variantevalall}{validation_report_header} =
      q?perl -nae' if ($_ =~/^ValidationReport\s+CompFeatureInput/ ) {print $_;last;}' ?;

    # Return ValidationReport and all line
    $regexp{variantevalall}{validation_report_data_all} =
q?perl -nae' if ( ($_ =~/^ValidationReport/) && ($_ =~/all\s/) && ($_ =~/none\s/)) {print $_;last;}' ?;

    # Return ValidationReport and known line
    $regexp{variantevalall}{validation_report_data_known} =
      q?perl -nae' if ( ($_ =~/^ValidationReport/) && ($_ =~/known\s/) ) {print $_;last;}' ?;

    # Return ValidationReport and novel line
    $regexp{variantevalall}{validation_report_data_novel} =
      q?perl -nae' if ( ($_ =~/^ValidationReport/) && ($_ =~/novel\s/) ) {print $_;last;}' ?;

    # Return VariantSummary and CompFeatureInput line from header
    $regexp{variantevalall}{variant_summary_header} =
      q?perl -nae' if ($_ =~/^VariantSummary\s+CompFeatureInput/ ) {print $_;last;}' ?;

    # Return VariantSummary and all line
    $regexp{variantevalall}{variant_summary_data_all} =
      q?perl -nae' if ( ($_ =~/^VariantSummary/) && ($_ =~/all\s/) ) {print $_;last;}' ?;

    # Return VariantSummary and known line
    $regexp{variantevalall}{variant_summary_data_known} =
      q?perl -nae' if ( ($_ =~/^VariantSummary/) && ($_ =~/known\s/) ) {print $_;last;}' ?;

    # Return VariantSummary and novel line
    $regexp{variantevalall}{variant_summary_data_novel} =
      q?perl -nae' if ( ($_ =~/^VariantSummary/) && ($_ =~/novel\s/) ) {print $_;last;}' ?;

    $regexp{variantevalexome} = $regexp{variantevalall};

    # Return varianteffectpredictor version
    $regexp{varianteffectpredictor}{version} =
      q?perl -nae 'if($_=~/##VEP="(\w+)"/) {print $1;last;}' ?;

    # Return varianteffectpredictor cache directory
    $regexp{varianteffectpredictor}{cache} =
      q?perl -nae 'if($_=~/##VEP=\w+\s+cache=(\S+)/) {print $1;last;}' ?;

    # Return varianteffectpredictor polyPhen version
    $regexp{varianteffectpredictor}{polyphen} =
      q?perl -nae 'if($_=~/##VEP=/ && $_=~/polyphen=(\S+)/) {print $1;last;}' ?;

    # Return varianteffectpredictor sift version
    $regexp{varianteffectpredictor}{sift} =
      q?perl -nae 'if($_=~/##VEP=/ && $_=~/sift=sift(\S+)/) {print $1;last;}' ?;

    # Return varianteffectpredictor geneBuild
    $regexp{varianteffectpredictor}{gene_build} =
      q?perl -nae 'if($_=~/##VEP=/ && $_=~/genebuild=(\S+)/) {print $1;last;}' ?;

    # Return varianteffectpredictor assembly
    $regexp{varianteffectpredictor}{assembly} =
      q?perl -nae 'if($_=~/##VEP=/ && $_=~/assembly=(\S+)/) {print $1;last;}' ?;

    # Return varianteffectpredictor HGMD-PUBLIC version
    $regexp{varianteffectpredictor}{hgmd_public} =
      q?perl -nae 'if($_=~/##VEP=/ && $_=~/HGMD-PUBLIC=(\S+)/) {print $1;last;}' ?;

    # Return varianteffectpredictor regbuild version
    $regexp{varianteffectpredictor}{reg_build} =
      q?perl -nae 'if($_=~/##VEP=/ && $_=~/regbuild=(\S+)/) {print $1;last;}' ?;

    # Return varianteffectpredictor gencode version
    $regexp{varianteffectpredictor}{gencode} =
      q?perl -nae 'if($_=~/##VEP=/ && $_=~/gencode=\S+\s+(\d+)/) {print $1;last;}' ?;

    # Return vcfparser version
    $regexp{vcfparser_ar}{version} =
      q?perl -nae 'if($_=~/##Software=<ID=mip,Version=(\d+.\d+.\d+)/) {print $1;last;}' ?;

    # Return sv_varianteffectpredictor version
    $regexp{sv_varianteffectpredictor}{version} =
      q?perl -nae 'if($_=~/##VEP="(\w+)"/) {print $1;last;}' ?;

    # Return sv_varianteffectpredictor cache directory
    $regexp{sv_varianteffectpredictor}{cache} =
      q?perl -nae 'if($_=~/##VEP=\w+\s+cache=(\S+)/) {print $1;last;}' ?;

    # Return sv_varianteffectpredictor polyPhen version
    $regexp{sv_varianteffectpredictor}{polyphen} =
      q?perl -nae 'if($_=~/##VEP=/ && $_=~/polyphen=(\S+)/) {print $1;last;}' ?;

    # Return sv_varianteffectpredictor sift version
    $regexp{sv_varianteffectpredictor}{sift} =
      q?perl -nae 'if($_=~/##VEP=/ && $_=~/sift=sift(\S+)/) {print $1;last;}' ?;

    # Return sv_varianteffectpredictor geneBuild
    $regexp{sv_varianteffectpredictor}{gene_build} =
      q?perl -nae 'if($_=~/##VEP=/ && $_=~/genebuild=(\S+)/) {print $1;last;}' ?;

    # Return sv_varianteffectpredictor assembly
    $regexp{sv_varianteffectpredictor}{assembly} =
      q?perl -nae 'if($_=~/##VEP=/ && $_=~/assembly=(\S+)/) {print $1;last;}' ?;

    # Return sv_varianteffectpredictor HGMD-PUBLIC version
    $regexp{sv_varianteffectpredictor}{hgmd_public} =
      q?perl -nae 'if($_=~/##VEP=/ && $_=~/HGMD-PUBLIC=(\S+)/) {print $1;last;}' ?;

    # Return sv_varianteffectpredictor regbuild version
    $regexp{sv_varianteffectpredictor}{reg_build} =
      q?perl -nae 'if($_=~/##VEP=/ && $_=~/regbuild=(\S+)/) {print $1;last;}' ?;

    # Return sv_varianteffectpredictor gencode version
    $regexp{sv_varianteffectpredictor}{gencode} =
      q?perl -nae 'if($_=~/##VEP=/ && $_=~/gencode=\S+\s+(\d+)/) {print $1;last;}' ?;

    # Return sv_vcfparser version
    $regexp{sv_vcfparser}{version} =
q?perl -nae 'if($_=~/##Software=<ID=mip,Version=(\d+.\d+.\d+)/) {print $1;last;} else { if($_=~/#CHROM/) {last;} }' ?;

    ## Writes a YAML hash to file
    write_to_file(
        {
            data_href => \%regexp,
            format    => q{yaml},
            path      => $print_regexp_outfile,
        }
    );
    $log->info( q{Wrote regexp YAML file to: } . $print_regexp_outfile );
    exit;
}

1;
