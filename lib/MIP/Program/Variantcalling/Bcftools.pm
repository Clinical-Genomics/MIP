package MIP::Program::Variantcalling::Bcftools;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use FindBin qw{ $Bin };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;

## CPANM
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.05;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ bcftools_annotate bcftools_call bcftools_concat bcftools_filter bcftools_index bcftools_merge bcftools_norm bcftools_reheader bcftools_rename_vcf_samples bcftools_roh bcftools_stats bcftools_view bcftools_view_and_index_vcf};

}

## Constants
Readonly my $COMMA        => q{,};
Readonly my $DOUBLE_QUOTE => q{"};
Readonly my $PIPE         => q{|};
Readonly my $NEWLINE      => qq{\n};
Readonly my $SPACE        => q{ };

sub bcftools_call {

## Function : Perl wrapper for writing bcftools call recipe to $FILEHANDLE or return commands array. Based on bcftools 1.3.1.
## Returns  : @commands
## Arguments: $form_fields_ref        => Output format fields {REF}
##          : $outfile_path           => Outfile path to write to
##          : $infile_path            => Infile path to read from
##          : $stderrfile_path        => Stderr file path to write to {OPTIONAL}
##          : $stderrfile_path_append => Append stderr info to file path
##          : $FILEHANDLE             => Filehandle to write to
##          : $samples_file           => PED file or a file with an optional column with sex
##          : $constrain              => One of: alleles, trio
##          : $multiallelic_caller    => Alternative model for multiallelic and rare-variant calling
##          : $output_type            => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##          : $variants_only          => Output variant sites only

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $form_fields_ref;
    my $outfile_path;
    my $infile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $FILEHANDLE;
    my $samples_file;
    my $constrain;

    ## Default(s)
    my $multiallelic_caller;
    my $output_type;
    my $variants_only;

    my $tmpl = {
        form_fields_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$form_fields_ref,
        },
        outfile_path    => { strict_type => 1, store => \$outfile_path, },
        infile_path     => { strict_type => 1, store => \$infile_path, },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path, },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append, },
        FILEHANDLE   => { store       => \$FILEHANDLE, },
        samples_file => { strict_type => 1, store => \$samples_file, },
        constrain => {
            allow       => [ undef, qw{ alleles trio } ],
            strict_type => 1,
            store       => \$constrain,
        },
        multiallelic_caller => {
            default     => 1,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$multiallelic_caller,
        },
        output_type => {
            default     => q{v},
            allow       => [qw{ b u z v }],
            strict_type => 1,
            store       => \$output_type,
        },
        variants_only => {
            default     => 1,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$variants_only,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ bcftools call };

    ## Options
    if ($multiallelic_caller) {

        push @commands, q{--multiallelic-caller};
    }

    if ( @{$form_fields_ref} ) {

        push @commands, q{--format-fields} . $SPACE . join $COMMA,
          @{$form_fields_ref};
    }

    if ($variants_only) {

        push @commands, q{--variants-only};
    }

    if ($samples_file) {

        push @commands, q{--samples-file} . $SPACE . $samples_file;
    }

    if ($constrain) {

        push @commands, q{--constrain} . $SPACE . $constrain;
    }

    if ($output_type) {

        #Specify output type
        push @commands, q{--output-type} . $SPACE . $output_type;
    }

    if ($outfile_path) {

        #Specify output filename
        push @commands, q{--output} . $SPACE . $outfile_path;
    }

    ## Infile
    if ($infile_path) {

        push @commands, $infile_path;
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;

}

sub bcftools_index {

## Function : Perl wrapper for writing bcftools index recipe to $FILEHANDLE or return commands array. Based on bcftools 1.3.1.
## Returns  : @commands
## Arguments: $infile_path            => Infile path to read from
##          : $stderrfile_path        => Stderr file path to write to {OPTIONAL}
##          : $stderrfile_path_append => Append stderr info to file path
##          : $FILEHANDLE             => Filehandle to write to
##          : $output_type            => 'csi' generate CSI-format index, 'tbi' generate TBI-format index

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $FILEHANDLE;

    ## Default(s)
    my $output_type;

    my $tmpl = {
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path,
        },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path, },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append, },
        FILEHANDLE  => { store => \$FILEHANDLE, },
        output_type => {
            default     => q{csi},
            allow       => [qw{ csi tbi }],
            strict_type => 1,
            store       => \$output_type,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ bcftools index };

    ## Options
    if ($output_type) {

        #Specify output type
        push @commands, q{--} . $output_type;
    }

    ## Infile
    if ($infile_path) {

        push @commands, $infile_path;
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;

}

sub bcftools_view {

## Function : Perl wrapper for writing bcftools view recipe to $FILEHANDLE or return commands array. Based on bcftools 1.4.1.
## Returns  : @commands
## Arguments: $apply_filters_ref      => Require at least one of the listed FILTER strings
##          : $exclude_types_ref      => Exclude comma-separated list of variant types: snps,indels,mnps,other
##          : $exclude                => Exclude sites for which the expression is true
##          : $include                => Include only sites for which the expression is true
##          : $sample                 => Comma separated list of samples to include (or exclude with "^" prefix)
##          : $infile_path            => Infile path to read from
##          : $outfile_path           => Outfile path to write to
##          : $stderrfile_path        => Stderr file path to write to
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile file path to write to
##          : $FILEHANDLE             => Filehandle to write to
##          : $output_type            => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $apply_filters_ref;
    my $exclude_types_ref;
    my $exclude;
    my $include;
    my $sample;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $FILEHANDLE;

    ## Default(s)
    my $output_type;

    my $tmpl = {
        apply_filters_ref =>
          { default => [], strict_type => 1, store => \$apply_filters_ref, },
        exclude_types_ref =>
          { default => [], strict_type => 1, store => \$exclude_types_ref, },
        exclude         => { strict_type => 1, store => \$exclude, },
        include         => { strict_type => 1, store => \$include, },
        sample          => { strict_type => 1, store => \$sample, },
        infile_path     => { strict_type => 1, store => \$infile_path, },
        outfile_path    => { strict_type => 1, store => \$outfile_path, },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path, },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append, },
        stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path, },
        FILEHANDLE  => { store => \$FILEHANDLE, },
        output_type => {
            default     => q{v},
            allow       => [qw{ b u z v }],
            strict_type => 1,
            store       => \$output_type,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ bcftools view };

    ## Options
    if ( @{$apply_filters_ref} ) {

        push @commands, q{--apply-filters} . $SPACE . join $COMMA,
          @{$apply_filters_ref};
    }

    if ( @{$exclude_types_ref} ) {

        push @commands, q{--exclude-types} . $SPACE . join $COMMA,
          @{$exclude_types_ref};
    }

    if ($exclude) {

        push @commands, q{--exclude} . $SPACE . $exclude;
    }

    if ($include) {

        push @commands, q{--include} . $SPACE . $include;
    }

    if ($sample) {

        push @commands, q{--samples} . $SPACE . $sample;
    }
    if ($output_type) {

        #Specify output type
        push @commands, q{--output-type} . $SPACE . $output_type;
    }

    if ($outfile_path) {

        #Specify output filename
        push @commands, q{--output-file} . $SPACE . $outfile_path;
    }

    ## Infile
    if ($infile_path) {

        push @commands, $infile_path;
    }

    push @commands,
      unix_standard_streams(
        {
            stdoutfile_path        => $stdoutfile_path,
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;

}

sub bcftools_filter {

## Function : Perl wrapper for writing bcftools filter recipe to $FILEHANDLE or return commands array. Based on bcftools 1.4.1.
## Returns  : @commands
## Arguments: $infile_path            => Infile paths
##          : $stdoutfile_path        => Stdoutfile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $FILEHANDLE             => Filehandle to write to
##          : $exclude                => Exclude sites for which the expression is true
##          : $include                => Include only sites for which the expression is true
##          : $soft_filter            => Annotate FILTER column with <string> or unique filter name
##          : $snp_gap                => Filter SNPs within <int> base pairs of an indel
##          : $indel_gap              => Filter clusters of indels separated by <int> or fewer base pairs allowing only one to pass

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $stdoutfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $FILEHANDLE;
    my $exclude;
    my $include;
    my $soft_filter;
    my $snp_gap;
    my $indel_gap;

    my $tmpl = {
        infile_path     => { strict_type => 1, store => \$infile_path, },
        stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path, },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path, },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append, },
        FILEHANDLE  => { store       => \$FILEHANDLE, },
        exclude     => { strict_type => 1, store => \$exclude, },
        include     => { strict_type => 1, store => \$include, },
        soft_filter => { strict_type => 1, store => \$soft_filter, },
        snp_gap => {
            allow       => qr/ ^\d+$ /sxm,
            strict_type => 1,
            store       => \$snp_gap,
        },
        indel_gap => {
            allow       => qr/ ^\d+$ /sxm,
            strict_type => 1,
            store       => \$indel_gap,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ bcftools filter };

    ## Options
    if ($exclude) {

        push @commands, q{--exclude} . $SPACE . $exclude;
    }
    if ($include) {

        push @commands, q{--include} . $SPACE . $include;
    }
    if ($soft_filter) {

        push @commands, q{--soft-filter} . $SPACE . $soft_filter;
    }

    if ($snp_gap) {

        push @commands, q{--SnpGap} . $SPACE . $snp_gap;
    }

    if ($indel_gap) {

        push @commands, q{--IndelGap} . $SPACE . $indel_gap;
    }

    ## Infile
    if ($infile_path) {

        push @commands, $infile_path;
    }

    push @commands,
      unix_standard_streams(
        {
            stdoutfile_path        => $stdoutfile_path,
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;

}

sub bcftools_norm {

## Function : Perl wrapper for writing bcftools norm recipe to $FILEHANDLE or return commands array. Based on bcftools 1.3.1.
## Returns  : @commands
## Arguments: $outfile_path           => Outfile path to write to
##          : $reference_path         => Human genome reference path
##          : $infile_path            => Infile path to read from
##          : $stderrfile_path        => Stderr file path to write to {OPTIONAL}
##          : $stderrfile_path_append => Append stderr info to file path
##          : $FILEHANDLE             => Filehandle to write to
##          : $multiallelic           => To split/join multiallelic calls or not
##          : $output_type            => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##          : $multiallelic_type      => Type of multiallelic to split/join {OPTIONAL}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $outfile_path;
    my $reference_path;
    my $infile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $FILEHANDLE;
    my $multiallelic;

    ## Default(s)
    my $output_type;
    my $multiallelic_type;

    my $tmpl = {
        outfile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfile_path,
        },
        reference_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$reference_path,
        },
        infile_path     => { strict_type => 1, store => \$infile_path, },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path, },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append, },
        FILEHANDLE   => { store => \$FILEHANDLE, },
        multiallelic => {
            allow       => [qw{ + - }],
            strict_type => 1,
            store       => \$multiallelic,
        },
        output_type => {
            default     => q{v},
            allow       => [qw{ b u z v }],
            strict_type => 1,
            store       => \$output_type,
        },
        multiallelic_type => {
            default     => q{both},
            allow       => [qw{ snps indels both any }],
            strict_type => 1,
            store       => \$multiallelic_type,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ bcftools norm };

    ## Options
    if ($multiallelic) {

        push @commands,
          q{--multiallelics} . $SPACE . $multiallelic . $multiallelic_type;
    }

    if ($reference_path) {

        push @commands, q{--fasta-ref} . $SPACE . $reference_path;
    }

    if ($output_type) {

        #Specify output type
        push @commands, q{--output-type} . $SPACE . $output_type;
    }

    if ($outfile_path) {

        #Specify output filename
        push @commands, q{--output} . $SPACE . $outfile_path;
    }

    ## Infile
    if ($infile_path) {

        push @commands, $infile_path;
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;

}

sub bcftools_merge {

## Function : Perl wrapper for writing bcftools merge recipe to $FILEHANDLE or return commands array. Based on bcftools 1.3.1.
## Returns  : @commands
## Arguments: $infile_paths_ref       => Infile path to read from
##          : $outfile_path           => Outfile path to write to
##          : $stderrfile_path        => Stderr file path to write to
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile file path to write to
##          : $FILEHANDLE             => Filehandle to write to
##          : $output_type            => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_paths_ref;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $FILEHANDLE;

    ## Default(s)
    my $output_type;

    my $tmpl = {
        infile_paths_ref =>
          { default => [], strict_type => 1, store => \$infile_paths_ref, },
        outfile_path    => { strict_type => 1, store => \$outfile_path, },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path, },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append, },
        stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path, },
        FILEHANDLE  => { store => \$FILEHANDLE, },
        output_type => {
            default     => q{v},
            allow       => [qw{ b u z v}],
            strict_type => 1,
            store       => \$output_type,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ bcftools merge };

    ## Options
    if ($output_type) {

        #Specify output type
        push @commands, q{--output-type} . $SPACE . $output_type;
    }

    if ($outfile_path) {

        #Specify output filename
        push @commands, q{--output} . $SPACE . $outfile_path;
    }

    if ( @{$infile_paths_ref} ) {

        push @commands, join $SPACE, @{$infile_paths_ref};
    }

    push @commands,
      unix_standard_streams(
        {
            stdoutfile_path        => $stdoutfile_path,
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;
}

sub bcftools_concat {

## Function : Perl wrapper for writing bcftools concat recipe to $FILEHANDLE or return commands array. Based on bcftools 1.4.1.
## Returns  : @commands
## Arguments: $infile_paths_ref       => Infile path to read from
##          : $outfile_path           => Outfile path to write to
##          : $stderrfile_path        => Stderr file path to write to
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile file path to write to
##          : $FILEHANDLE             => Filehandle to write to
##          : $allow_overlaps         => First coordinate of the next file can precede last record of the current file
##          : $rm_dups      => Output duplicate records present in multiple files only once: <snps|indels|both|all|none>
    ##          : $output_type            => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_paths_ref;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $FILEHANDLE;

    ## Default(s)
    my $allow_overlaps;
    my $rm_dups;
    my $output_type;

    my $tmpl = {
        infile_paths_ref =>
          { default => [], strict_type => 1, store => \$infile_paths_ref, },
        outfile_path    => { strict_type => 1, store => \$outfile_path, },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path, },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append, },
        stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path, },
        FILEHANDLE => { store => \$FILEHANDLE, },
        rm_dups    => {
            default     => q{all},
            allow       => [qw{ snps indels both all none }],
            strict_type => 1,
            store       => \$rm_dups,
        },
        allow_overlaps => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$allow_overlaps,
        },
        output_type => {
            default     => q{v},
            allow       => [qw{ b u z v }],
            strict_type => 1,
            store       => \$output_type,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ bcftools concat };

    ## Options
    if ($allow_overlaps) {

        push @commands, q{--allow-overlaps};
    }

    if ($rm_dups) {

        push @commands, q{--rm-dups} . $SPACE . $rm_dups;
    }

    if ($output_type) {

        #Specify output type
        push @commands, q{--output-type} . $SPACE . $output_type;
    }

    if ($outfile_path) {

        #Specify output filename
        push @commands, q{--output} . $SPACE . $outfile_path;
    }

    ## Infile
    if ( @{$infile_paths_ref} ) {

        push @commands, join $SPACE, @{$infile_paths_ref};
    }

    push @commands,
      unix_standard_streams(
        {
            stdoutfile_path        => $stdoutfile_path,
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    return @commands;

}

sub bcftools_annotate {

## Function : Perl wrapper for writing bcftools annotate recipe to $FILEHANDLE or return commands array. Based on bcftools 1.3.1.
## Returns  : @commands
## Arguments: $remove_ids_ref         => List of annotations to remove
##          : $infile_path            => Infile path to read from
##          : $outfile_path           => Outfile path to write to
##          : $stderrfile_path        => Stderr file path to write to {OPTIONAL}
##          : $stderrfile_path_append => Append stderr info to file path
##          : $FILEHANDLE             => Filehandle to write to
##          : $samples_file           => File of samples to annotate
##          : $headerfile_path        => File with lines which should be appended to the VCF header
##          : $set_id                 => Set ID column
##          : $output_type            => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $remove_ids_ref;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $FILEHANDLE;
    my $samples_file;
    my $headerfile_path;
    my $set_id;

    ## Default(s)
    my $output_type;

    my $tmpl = {
        remove_ids_ref => {
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$remove_ids_ref,
        },
        infile_path     => { strict_type => 1, store => \$infile_path, },
        outfile_path    => { strict_type => 1, store => \$outfile_path, },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path, },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append, },
        FILEHANDLE      => { store       => \$FILEHANDLE, },
        samples_file    => { strict_type => 1, store => \$samples_file, },
        headerfile_path => { strict_type => 1, store => \$headerfile_path, },
        set_id          => { strict_type => 1, store => \$set_id, },
        output_type => {
            default     => q{v},
            allow       => [qw{ b u z v}],
            strict_type => 1,
            store       => \$output_type,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ bcftools annotate };

    ## Options
    if ( @{$remove_ids_ref} ) {

        push @commands, q{--remove} . $SPACE . join $COMMA, @{$remove_ids_ref};
    }

    if ($set_id) {

        push @commands, q{--set-id} . $SPACE . $set_id;
    }

    if ($samples_file) {

        push @commands, q{--samples-file} . $SPACE . $samples_file;
    }

    if ($headerfile_path) {

        push @commands, q{--header-lines} . $SPACE . $headerfile_path;
    }

    if ($output_type) {

        #Specify output type
        push @commands, q{--output-type} . $SPACE . $output_type;
    }

    ## Infile
    if ($infile_path) {

        push @commands, $infile_path;
    }

    if ($outfile_path) {

        #Specify output filename
        push @commands, q{>} . $SPACE . $outfile_path;
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;
}

sub bcftools_roh {

## Function : Perl wrapper for writing bcftools roh recipe to $FILEHANDLE or return commands array. Based on bcftools 1.4.1.
## Returns  : @commands
## Arguments: $sample_ids_ref         => Sample to analyze
##          : $infile_path            => Infile path to read from
##          : $outfile_path           => Outfile path to write to
##          : $stderrfile_path        => Stderr file path to write to
##          : $stdoutfile_path        => Stdoutfile file path to write to
##          : $stderrfile_path_append => Append stderr info to file path
##          : $FILEHANDLE             => Filehandle to write to
##          : $af_file_path           => Read allele frequencies from file (CHR\tPOS\tREF,ALT\tAF)
##          : $skip_indels            => Skip indels as their genotypes are enriched for errors

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_ids_ref;
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $FILEHANDLE;
    my $af_file_path;

    ## Default(s)
    my $skip_indels;

    my $tmpl = {
        sample_ids_ref =>
          { default => [], strict_type => 1, store => \$sample_ids_ref },
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path,
        },
        outfile_path    => { strict_type => 1, store => \$outfile_path, },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path, },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append, },
        stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path, },
        FILEHANDLE   => { store       => \$FILEHANDLE, },
        af_file_path => { strict_type => 1, store => \$af_file_path, },
        skip_indels => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$skip_indels,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ bcftools roh };

    ## Options
    if ($af_file_path) {

        push @commands, q{--AF-file} . $SPACE . $af_file_path;
    }

    if ($skip_indels) {

        push @commands, q{--skip-indels};
    }

    if ( @{$sample_ids_ref} ) {

        push @commands, q{--samples} . $SPACE . join $COMMA, @{$sample_ids_ref};
    }

    ## Infile
    if ($infile_path) {

        push @commands, $infile_path;
    }

    if ($outfile_path) {

        #Specify output filename
        push @commands, q{>} . $SPACE . $outfile_path;
    }

    push @commands,
      unix_standard_streams(
        {
            stdoutfile_path        => $stdoutfile_path,
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;

}

sub bcftools_stats {

## Function : Perl wrapper for writing bcftools stats recipe to already open $FILEHANDLE or return commands array. Based on bcftools 1.3.1.
## Returns  : @commands
## Arguments: $infile_path            => Infile path
##          : $outfile_path           => Outfile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $FILEHANDLE             => Filehandle to write to
##          : $append_stderr_info     => Append stderr info to file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $FILEHANDLE;

    ## Default(s)
    my $append_stderr_info;

    my $tmpl = {
        infile_path     => { strict_type => 1, store => \$infile_path, },
        outfile_path    => { strict_type => 1, store => \$outfile_path, },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path, },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append, },
        stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path, },
        FILEHANDLE         => { store => \$FILEHANDLE, },
        append_stderr_info => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$append_stderr_info,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ bcftools stats };

    ## Infile
    if ($infile_path) {

        push @commands, $infile_path;
    }

    ## Outfile
    if ($outfile_path) {

        push @commands, q{>} . $SPACE . $outfile_path;
    }

    if ( $append_stderr_info && $stderrfile_path ) {

        push @commands, q{2>>} . $SPACE . $stderrfile_path;
    }

    push @commands,
      unix_standard_streams(
        {
            stdoutfile_path        => $stdoutfile_path,
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;

}

sub bcftools_reheader {

## Function : Perl wrapper for writing bcftools reheader recipe to already open $FILEHANDLE or return commands array. Based on bcftools 1.3.1.
## Returns  : @commands
## Arguments: $infile_path            => Infile path
##          : $outfile_path           => Outfile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $samples_file_path      => File with new sample names
##          : $FILEHANDLE             => Filehandle to write to
##          : $append_stderr_info     => Append stderr info to file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $outfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $samples_file_path;
    my $FILEHANDLE;

    ## Default(s)
    my $append_stderr_info;

    my $tmpl = {
        infile_path     => { strict_type => 1, store => \$infile_path, },
        outfile_path    => { strict_type => 1, store => \$outfile_path, },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path, },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append, },
        stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path, },
        samples_file_path =>
          { strict_type => 1, store => \$samples_file_path, },
        FILEHANDLE         => { store => \$FILEHANDLE, },
        append_stderr_info => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$append_stderr_info,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ bcftools reheader };

    ## Options
    if ($samples_file_path) {

        push @commands, q{--samples} . $SPACE . $samples_file_path;
    }

    ## Infile
    if ($infile_path) {

        push @commands, $infile_path;
    }

    ## Outfile
    if ($outfile_path) {

        push @commands, q{>} . $SPACE . $outfile_path;
    }

    if ( $append_stderr_info && $stderrfile_path ) {

        push @commands, q{2>>} . $SPACE . $stderrfile_path;
    }

    push @commands,
      unix_standard_streams(
        {
            stdoutfile_path        => $stdoutfile_path,
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            separator    => $SPACE,
            FILEHANDLE   => $FILEHANDLE,
        }
    );

    return @commands;

}

sub bcftools_rename_vcf_samples {

## Function : Rename vcf samples. The samples array will replace the sample names in the same order as supplied.
## Returns  :
## Arguments: $sample_ids_ref => Samples to rename in the same order as in the vcf {REF}
##          : $temp_directory => Temporary directory
##          : $infile         => Vcf infile to rename samples for
##          : $outfile        => Output vcf with samples renamed
##          : $FILEHANDLE     => Filehandle to write to
##          : $output_type    => Output type

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sample_ids_ref;
    my $temp_directory;
    my $infile;
    my $outfile;
    my $FILEHANDLE;

    ## Default(s)
    my $output_type;

    my $tmpl = {
        sample_ids_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$sample_ids_ref,
        },
        temp_directory => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$temp_directory,
        },
        infile =>
          { required => 1, defined => 1, strict_type => 1, store => \$infile, },
        outfile => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfile,
        },
        FILEHANDLE  => { required => 1, defined => 1, store => \$FILEHANDLE, },
        output_type => {
            default     => q{v},
            allow       => [qw{ b u z v }],
            strict_type => 1,
            store       => \$output_type,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Gnu::Coreutils qw{ gnu_printf };

    ## Create new sample names file
    say {$FILEHANDLE} q{## Create new sample(s) names file};

    ## Get parameters
    my $format_string = $DOUBLE_QUOTE;
  SAMPLE_ID:
    foreach my $sample_id ( @{$sample_ids_ref} ) {

        $format_string .= $sample_id . q{\n};
    }
    $format_string .= $DOUBLE_QUOTE;
    gnu_printf(
        {
            format_string   => $format_string,
            stdoutfile_path => catfile( $temp_directory, q{sample_name.txt} ),
            FILEHANDLE      => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    ## Rename samples in VCF
    say {$FILEHANDLE} q{## Rename sample(s) names in VCF file};
    bcftools_reheader(
        {
            infile_path       => $infile,
            samples_file_path => catfile( $temp_directory, q{sample_name.txt} ),
            FILEHANDLE        => $FILEHANDLE,
        }
    );
    ## Pipe
    print {$FILEHANDLE} $PIPE . $SPACE;

    bcftools_view(
        {
            outfile_path => $outfile,
            output_type  => q{v},
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;
    return;
}

sub bcftools_view_and_index_vcf {

## Function : View variant calling file and index.
## Returns  :
## Arguments: $infile_path         => Path to infile to compress and index
##          : $FILEHANDLE          => SBATCH script FILEHANDLE to print to
##          : $outfile_path_prefix => Out file path no file_ending {Optional}
##          : $output_type         => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##          : $index               => Generate index of reformated file
##          : $index_type          => Type of index

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $FILEHANDLE;
    my $outfile_path_prefix;

    ## Default(s)
    my $output_type;
    my $index;
    my $index_type;

    my $tmpl = {
        infile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$infile_path,
        },
        FILEHANDLE => { required => 1, defined => 1, store => \$FILEHANDLE, },
        outfile_path_prefix =>
          { strict_type => 1, store => \$outfile_path_prefix, },
        output_type => {
            default     => q{v},
            allow       => [qw{ b u z v }],
            strict_type => 1,
            store       => \$output_type,
        },
        index => {
            default     => 1,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$index,
        },
        index_type => {
            default     => q{csi},
            allow       => [ undef, qw{ csi tbi } ],
            strict_type => 1,
            store       => \$index_type,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $outfile_path;
    my %output_type_ending = (
        b => q{.bcf},
        u => q{.bcf},
        z => q{.vcf.gz},
        v => q{.vcf},
    );

    if ( defined $outfile_path_prefix ) {

        $outfile_path =
          $outfile_path_prefix . $output_type_ending{$output_type};
    }

    say {$FILEHANDLE} q{## Reformat variant calling file};

    bcftools_view(
        {
            infile_path  => $infile_path,
            outfile_path => $outfile_path,
            output_type  => $output_type,
            FILEHANDLE   => $FILEHANDLE,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    if ($index) {

        say {$FILEHANDLE} q{## Index};

        bcftools_index(
            {
                infile_path => $outfile_path,
                output_type => $index_type,
                FILEHANDLE  => $FILEHANDLE,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }
    return;
}

1;
