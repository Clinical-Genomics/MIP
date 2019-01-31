package MIP::Program::Variantcalling::Bcftools;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Program::Base::Bcftools qw{ bcftools_base };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.13;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      bcftools_annotate
      bcftools_call
      bcftools_concat
      bcftools_create_reheader_samples_file
      bcftools_filter
      bcftools_index
      bcftools_merge
      bcftools_mpileup
      bcftools_norm
      bcftools_reheader
      bcftools_rename_vcf_samples
      bcftools_roh
      bcftools_stats
      bcftools_view
      bcftools_view_and_index_vcf
    };
}

## Constants
Readonly my $COMMA        => q{,};
Readonly my $DOUBLE_QUOTE => q{"};
Readonly my $PIPE         => q{|};
Readonly my $NEWLINE      => qq{\n};
Readonly my $SPACE        => q{ };

sub bcftools_annotate {

## Function : Perl wrapper for writing bcftools annotate recipe to $FILEHANDLE or return commands array. Based on bcftools 1.9.
## Returns  : @commands
## Arguments: $annotations_file_path  => VCF file or tabix-indexed file path with annotations: CHR\tPOS[\tVALUE]+
##          : $columns_name           => List of columns in the annotation file, e.g. CHROM,POS,REF,ALT,-,INFO/TAG
##          : $FILEHANDLE             => Filehandle to write to
##          : $headerfile_path        => File with lines which should be appended to the VCF header
##          : $infile_path            => Infile path to read from
##          : $outfile_path           => Outfile path to write to
##          : $output_type            => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##          : $regions_ref            => Regions to process {REF}
##          : $remove_ids_ref         => List of annotations to remove
##          : $samples_file_path      => File of samples to annotate
##          : $samples_ref            => Samples to include or exclude if prefixed with "^"
##          : $set_id                 => Set ID column
##          : $stderrfile_path        => Stderr file path to write to {OPTIONAL}
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $annotations_file_path;
    my $columns_name;
    my $FILEHANDLE;
    my $infile_path;
    my $headerfile_path;
    my $outfile_path;
    my $regions_ref;
    my $remove_ids_ref;
    my $samples_file_path;
    my $samples_ref;
    my $set_id;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $output_type;

    my $tmpl = {
        annotations_file_path => { store => \$annotations_file_path, strict_type => 1, },
        columns_name          => { store => \$columns_name,          strict_type => 1, },
        FILEHANDLE            => { store => \$FILEHANDLE, },
        headerfile_path       => { store => \$headerfile_path,       strict_type => 1, },
        infile_path           => { store => \$infile_path,           strict_type => 1, },
        outfile_path          => { store => \$outfile_path,          strict_type => 1, },
        output_type           => {
            allow       => [qw{ b u z v}],
            default     => q{v},
            store       => \$output_type,
            strict_type => 1,
        },
        regions_ref    => { default => [], store => \$regions_ref, strict_type => 1, },
        remove_ids_ref => {
            default     => [],
            defined     => 1,
            store       => \$remove_ids_ref,
            strict_type => 1,
        },
        samples_file_path => { store => \$samples_file_path, strict_type => 1, },
        samples_ref       => {
            default     => [],
            store       => \$samples_ref,
            strict_type => 1,
        },
        set_id          => { store => \$set_id,          strict_type => 1, },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ bcftools annotate };

    ## Bcftools base args
    @commands = bcftools_base(
        {
            commands_ref      => \@commands,
            outfile_path      => $outfile_path,
            output_type       => $output_type,
            regions_ref       => $regions_ref,
            samples_file_path => $samples_file_path,
            samples_ref       => $samples_ref,
        }
    );

    ## Options
    if ($annotations_file_path) {

        push @commands, q{--annotations} . $SPACE . $annotations_file_path;
    }

    if ($columns_name) {

        push @commands, q{--columns} . $SPACE . $columns_name;
    }

    if ( @{$remove_ids_ref} ) {

        push @commands, q{--remove} . $SPACE . join $COMMA, @{$remove_ids_ref};
    }

    if ($set_id) {

        push @commands, q{--set-id} . $SPACE . $set_id;
    }

    if ($headerfile_path) {

        push @commands, q{--header-lines} . $SPACE . $headerfile_path;
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
            stdoutfile_path        => $stdoutfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );

    return @commands;
}

sub bcftools_call {

## Function : Perl wrapper for writing bcftools call recipe to $FILEHANDLE or return commands array. Based on bcftools 1.6.
## Returns  : @commands
## Arguments: $constrain              => One of: alleles, trio
##          : $FILEHANDLE             => Filehandle to write to
##          : $form_fields_ref        => Output format fields {REF}
##          : $infile_path            => Infile path to read from
##          : $multiallelic_caller    => Alternative model for multiallelic and rare-variant calling
##          : $outfile_path           => Outfile path to write to
##          : $output_type            => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##          : $regions_ref            => Regions to process {REF}
##          : $samples_file_path      => PED file or a file with an optional column with sex
##          : $samples_ref            => Samples to include or exclude if prefixed with "^"
##          : $stderrfile_path        => Stderr file path to write to {OPTIONAL}
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $variants_only          => Output variant sites only

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $constrain;
    my $FILEHANDLE;
    my $form_fields_ref;
    my $infile_path;
    my $outfile_path;
    my $regions_ref;
    my $samples_file_path;
    my $samples_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $multiallelic_caller;
    my $output_type;
    my $variants_only;

    my $tmpl = {
        constrain => {
            allow       => [ undef, qw{ alleles trio } ],
            store       => \$constrain,
            strict_type => 1,
        },
        FILEHANDLE      => { store => \$FILEHANDLE, },
        form_fields_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$form_fields_ref,
            strict_type => 1,
        },
        infile_path         => { store => \$infile_path, strict_type => 1, },
        multiallelic_caller => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$multiallelic_caller,
            strict_type => 1,
        },
        outfile_path => { store => \$outfile_path, strict_type => 1, },
        output_type  => {
            allow       => [qw{ b u z v }],
            default     => q{v},
            store       => \$output_type,
            strict_type => 1,
        },
        regions_ref => { default => [], store => \$regions_ref, strict_type => 1, },
        samples_file_path => { store => \$samples_file_path, strict_type => 1, },
        samples_ref       => {
            default     => [],
            store       => \$samples_ref,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
        variants_only   => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$variants_only,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ bcftools call };

    ## Bcftools base args
    @commands = bcftools_base(
        {
            commands_ref      => \@commands,
            outfile_path      => $outfile_path,
            output_type       => $output_type,
            regions_ref       => $regions_ref,
            samples_file_path => $samples_file_path,
            samples_ref       => $samples_ref,
        }
    );

    ## Options
    if ($multiallelic_caller) {

        push @commands, q{--multiallelic-caller};
    }

    if ( @{$form_fields_ref} ) {

        push @commands, q{--format-fields} . $SPACE . join $COMMA, @{$form_fields_ref};
    }

    if ($variants_only) {

        push @commands, q{--variants-only};
    }

    if ($constrain) {

        push @commands, q{--constrain} . $SPACE . $constrain;
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
            stdoutfile_path        => $stdoutfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );

    return @commands;
}

sub bcftools_concat {

## Function : Perl wrapper for writing bcftools concat recipe to $FILEHANDLE or return commands array. Based on bcftools 1.6.
## Returns  : @commands
## Arguments: $allow_overlaps         => First coordinate of the next file can precede last record of the current file
##          : $FILEHANDLE             => Filehandle to write to
##          : $infile_paths_ref       => Infile path to read from
##          : $outfile_path           => Outfile path to write to
##          : $output_type            => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##          : $regions_ref            => Regions to process {REF}
##          : $rm_dups                => Output duplicate records present in multiple files only once: <snps|indels|both|all|none>
##          : $samples_file_path      => File of samples to annotate
##          : $samples_ref            => Samples to include or exclude if prefixed with "^"
##          : $stderrfile_path        => Stderr file path to write to
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile file path to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_paths_ref;
    my $outfile_path;
    my $regions_ref;
    my $samples_file_path;
    my $samples_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $allow_overlaps;
    my $output_type;
    my $rm_dups;

    my $tmpl = {
        allow_overlaps => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$allow_overlaps,
            strict_type => 1,
        },
        FILEHANDLE => { store => \$FILEHANDLE, },
        infile_paths_ref =>
          { default => [], store => \$infile_paths_ref, strict_type => 1, },
        outfile_path => { store => \$outfile_path, strict_type => 1, },
        output_type  => {
            allow       => [qw{ b u z v }],
            default     => q{v},
            store       => \$output_type,
            strict_type => 1,
        },
        regions_ref => { default => [], store => \$regions_ref, strict_type => 1, },
        rm_dups     => {
            allow       => [qw{ 0 snps indels both all none }],
            default     => q{all},
            store       => \$rm_dups,
            strict_type => 1,
        },
        samples_file_path => { store => \$samples_file_path, strict_type => 1, },
        samples_ref       => {
            default     => [],
            store       => \$samples_ref,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ bcftools concat };

    ## Bcftools base args
    @commands = bcftools_base(
        {
            commands_ref      => \@commands,
            outfile_path      => $outfile_path,
            output_type       => $output_type,
            regions_ref       => $regions_ref,
            samples_file_path => $samples_file_path,
            samples_ref       => $samples_ref,
        }
    );

    ## Options
    if ($allow_overlaps) {

        push @commands, q{--allow-overlaps};
    }

    if ($rm_dups) {

        push @commands, q{--rm-dups} . $SPACE . $rm_dups;
    }

    ## Infile
    if ( @{$infile_paths_ref} ) {

        push @commands, join $SPACE, @{$infile_paths_ref};
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
            stdoutfile_path        => $stdoutfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );
    return @commands;

}

sub bcftools_filter {

## Function : Perl wrapper for writing bcftools filter recipe to $FILEHANDLE or return commands array. Based on bcftools 1.6.
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $exclude                => Exclude sites for which the expression is true
##          : $filter_mode            => "+": do not replace but add to existing FILTER; "x": reset filters at sites which pass
##          : $include                => Include only sites for which the expression is true
##          : $indel_gap              => Filter clusters of indels separated by <int> or fewer base pairs allowing only one to pass
##          : $infile_path            => Infile paths
##          : $outfile_path           => Outfile path to write to
##          : $output_type            => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##          : $regions_ref            => Regions to process {REF}
##          : $samples_file_path      => File of samples to annotate
##          : $samples_ref            => Samples to include or exclude if prefixed with "^"
##          : $soft_filter            => Annotate FILTER column with <string> or unique filter name
##          : $snp_gap                => Filter SNPs within <int> base pairs of an indel
##          : $stdoutfile_path        => Stdoutfile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $exclude;
    my $filter_mode;
    my $FILEHANDLE;
    my $include;
    my $indel_gap;
    my $infile_path;
    my $outfile_path;
    my $regions_ref;
    my $samples_file_path;
    my $samples_ref;
    my $soft_filter;
    my $snp_gap;
    my $stdoutfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;

    ## Default(s)
    my $output_type;

    my $tmpl = {
        FILEHANDLE => { store => \$FILEHANDLE, },
        exclude    => {
            store       => \$exclude,
            strict_type => 1,
        },
        filter_mode => {
            allow       => [qw{ + x }],
            default     => q{+},
            store       => \$filter_mode,
            strict_type => 1,
        },
        infile_path => {
            store       => \$infile_path,
            strict_type => 1,
        },
        include => {
            store       => \$include,
            strict_type => 1,
        },
        indel_gap => {
            allow       => qr/ ^\d+$ /sxm,
            store       => \$indel_gap,
            strict_type => 1,
        },
        outfile_path => {
            store       => \$outfile_path,
            strict_type => 1,
        },
        output_type => {
            allow       => [qw{ b u z v}],
            default     => q{v},
            store       => \$output_type,
            strict_type => 1,
        },
        regions_ref => {
            default     => [],
            store       => \$regions_ref,
            strict_type => 1,
        },
        samples_file_path => {
            store       => \$samples_file_path,
            strict_type => 1,
        },
        samples_ref => {
            default     => [],
            store       => \$samples_ref,
            strict_type => 1,
        },
        soft_filter => {
            store       => \$soft_filter,
            strict_type => 1,
        },
        snp_gap => {
            allow       => qr/ ^\d+$ /sxm,
            store       => \$snp_gap,
            strict_type => 1,
        },
        stderrfile_path => {
            store       => \$stderrfile_path,
            strict_type => 1,
        },
        stderrfile_path_append => {
            store       => \$stderrfile_path_append,
            strict_type => 1,
        },
        stdoutfile_path => {
            store       => \$stdoutfile_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ bcftools filter };

    ## Bcftools base args
    @commands = bcftools_base(
        {
            commands_ref      => \@commands,
            outfile_path      => $outfile_path,
            output_type       => $output_type,
            regions_ref       => $regions_ref,
            samples_file_path => $samples_file_path,
            samples_ref       => $samples_ref,
        }
    );

    ## Options
    if ($exclude) {

        push @commands, q{--exclude} . $SPACE . $exclude;
    }

    if ($filter_mode) {

        push @commands, q{--mode} . $SPACE . $filter_mode;
    }

    if ($include) {

        push @commands, q{--include} . $SPACE . $include;
    }

    if ($indel_gap) {

        push @commands, q{--IndelGap} . $SPACE . $indel_gap;
    }

    ## Infile
    if ($infile_path) {

        push @commands, $infile_path;
    }

    if ($soft_filter) {

        push @commands, q{--soft-filter} . $SPACE . $soft_filter;
    }

    if ($snp_gap) {

        push @commands, q{--SnpGap} . $SPACE . $snp_gap;
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
            stdoutfile_path        => $stdoutfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );

    return @commands;

}

sub bcftools_index {

## Function : Perl wrapper for writing bcftools index recipe to $FILEHANDLE or return commands array. Based on bcftools 1.6.
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $infile_path            => Infile path to read from
##          : $outfile_path           => Outfile path to write to
##          : $output_type            => 'csi' generate CSI-format index, 'tbi' generate TBI-format index
##          : $regions_ref            => Regions to process {REF}
##          : $samples_file_path      => File of samples to annotate
##          : $samples_ref            => Samples to include or exclude if prefixed with "^"
##          : $stderrfile_path        => Stderr file path to write to {OPTIONAL}
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_path;
    my $regions_ref;
    my $samples_file_path;
    my $samples_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $output_type;

    my $tmpl = {
        FILEHANDLE  => { store => \$FILEHANDLE, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outfile_path => { store => \$outfile_path, strict_type => 1, },
        output_type  => {
            allow       => [qw{ csi tbi }],
            default     => q{csi},
            store       => \$output_type,
            strict_type => 1,
        },
        regions_ref => { default => [], store => \$regions_ref, strict_type => 1, },
        samples_file_path => { store => \$samples_file_path, strict_type => 1, },
        samples_ref       => {
            default     => [],
            store       => \$samples_ref,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ bcftools index };

    ## Bcftools base args
    @commands = bcftools_base(
        {
            commands_ref      => \@commands,
            outfile_path      => $outfile_path,
            regions_ref       => $regions_ref,
            samples_file_path => $samples_file_path,
            samples_ref       => $samples_ref,
        }
    );

    ## Options
    # Special case: 'csi' or 'tbi'
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
            stdoutfile_path        => $stdoutfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );

    return @commands;

}

sub bcftools_merge {

## Function : Perl wrapper for writing bcftools merge recipe to $FILEHANDLE or return commands array. Based on bcftools 1.6.
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $infile_paths_ref       => Infile path to read from
##          : $outfile_path           => Outfile path to write to
##          : $output_type            => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##          : $regions_ref            => Regions to process {REF}
##          : $samples_file_path      => File of samples to annotate
##          : $samples_ref            => Samples to include or exclude if prefixed with "^"
##          : $stderrfile_path        => Stderr file path to write to
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile file path to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_paths_ref;
    my $outfile_path;
    my $regions_ref;
    my $samples_file_path;
    my $samples_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $output_type;

    my $tmpl = {
        FILEHANDLE => { store => \$FILEHANDLE, },
        infile_paths_ref =>
          { default => [], store => \$infile_paths_ref, strict_type => 1, },
        outfile_path => { store => \$outfile_path, strict_type => 1, },
        output_type  => {
            allow       => [qw{ b u z v}],
            default     => q{v},
            store       => \$output_type,
            strict_type => 1,
        },
        regions_ref => { default => [], store => \$regions_ref, strict_type => 1, },
        samples_file_path => { store => \$samples_file_path, strict_type => 1, },
        samples_ref       => {
            default     => [],
            store       => \$samples_ref,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ bcftools merge };

    ## Bcftools base args
    @commands = bcftools_base(
        {
            commands_ref      => \@commands,
            outfile_path      => $outfile_path,
            output_type       => $output_type,
            regions_ref       => $regions_ref,
            samples_file_path => $samples_file_path,
            samples_ref       => $samples_ref,
        }
    );

    ## Options
    if ( @{$infile_paths_ref} ) {

        push @commands, join $SPACE, @{$infile_paths_ref};
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
            stdoutfile_path        => $stdoutfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );

    return @commands;
}

sub bcftools_mpileup {

## Function : Perl wrapper for writing bcftools mpileup recipe to $FILEHANDLE. Based on bcftools 1.6 (using htslib 1.6).
## Returns  : @commands
##          : $adjust_mq                        => Adjust mapping quality
##          : $FILEHANDLE                       => Sbatch filehandle to write to
##          : $infile_paths_ref                 => Infile paths {REF}
##          : $outfile_path                     => Outfile path
##          : $output_tags_ref                  => Optional tags to output {REF}
##          : $output_type                      => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##          : $per_sample_increased_sensitivity => Apply -m and -F per-sample for increased sensitivity
##          : $referencefile_path               => Reference sequence file
##          : $regions_ref                      => Regions to process {REF}
##          : $samples_file_path                => File of samples to annotate
##          : $samples_ref                      => Samples to include or exclude if prefixed with "^"
##          : $stderrfile_path                  => Stderrfile path
##          : $stderrfile_path_append           => Stderrfile path append
##          : $stdoutfile_path                  => Stdoutfile file path to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_paths_ref;
    my $outfile_path;
    my $output_tags_ref;
    my $referencefile_path;
    my $regions_ref;
    my $samples_file_path;
    my $samples_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $adjust_mq;
    my $per_sample_increased_sensitivity;
    my $output_type;

    ## Constants
    Readonly my $ADJUST_MAPPING_QUALITY => 50;

    my $tmpl = {
        adjust_mq => {
            allow       => qr/ ^\d+$ /sxm,
            default     => $ADJUST_MAPPING_QUALITY,
            store       => \$adjust_mq,
            strict_type => 1,
        },
        FILEHANDLE       => { store => \$FILEHANDLE, },
        infile_paths_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$infile_paths_ref,
            strict_type => 1,
        },
        outfile_path    => { store => \$outfile_path, strict_type => 1, },
        output_tags_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$output_tags_ref,
            strict_type => 1,
        },
        output_type => {
            allow       => [qw{ b u z v}],
            default     => q{b},
            store       => \$output_type,
            strict_type => 1,
        },
        per_sample_increased_sensitivity => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$per_sample_increased_sensitivity,
            strict_type => 1,
        },
        referencefile_path => {
            defined     => 1,
            required    => 1,
            store       => \$referencefile_path,
            strict_type => 1,
        },
        regions_ref => { default => [], store => \$regions_ref, strict_type => 1, },
        samples_file_path => { store => \$samples_file_path, strict_type => 1, },
        samples_ref       => {
            default     => [],
            store       => \$samples_ref,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Array @commands stores commands depending on input parameters
    my @commands = qw{ bcftools mpileup };

    ## Bcftools base args
    @commands = bcftools_base(
        {
            commands_ref      => \@commands,
            regions_ref       => $regions_ref,
            outfile_path      => $outfile_path,
            output_type       => $output_type,
            samples_file_path => $samples_file_path,
            samples_ref       => $samples_ref,
        }
    );

    ## Options
    push @commands, q{--adjust-MQ} . $SPACE . $adjust_mq;

    if ($per_sample_increased_sensitivity) {

        push @commands, q{--per-sample-mF};
    }

    if ( @{$output_tags_ref} ) {

        push @commands, q{--annotate} . $SPACE . join $COMMA, @{$output_tags_ref};
    }

    # Reference sequence file
    push @commands, q{--fasta-ref} . $SPACE . $referencefile_path;

    ## Infile
    push @commands, join $SPACE, @{$infile_paths_ref};

    # Redirect stderr output to program specific stderr file
    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
            stdoutfile_path        => $stdoutfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub bcftools_norm {

## Function : Perl wrapper for writing bcftools norm recipe to $FILEHANDLE or return commands array. Based on bcftools 1.6.
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $infile_path            => Infile path to read from
##          : $multiallelic           => To split/join multiallelic calls or not
##          : $multiallelic_type      => Type of multiallelic to split/join {OPTIONAL}
##          : $outfile_path           => Outfile path to write to
##          : $output_type            => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##          : $reference_path         => Human genome reference path
##          : $regions_ref            => Regions to process {REF}
##          : $samples_file_path      => File of samples to annotate
##          : $samples_ref            => Samples to include or exclude if prefixed with "^"
##          : $stderrfile_path        => Stderr file path to write to {OPTIONAL}
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $multiallelic;
    my $outfile_path;
    my $reference_path;
    my $regions_ref;
    my $samples_file_path;
    my $samples_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $multiallelic_type;
    my $output_type;

    my $tmpl = {
        FILEHANDLE  => { store => \$FILEHANDLE, },
        infile_path => { store => \$infile_path, strict_type => 1, },
        multiallelic => {
            allow       => [qw{ + - }],
            store       => \$multiallelic,
            strict_type => 1,
        },
        multiallelic_type => {
            allow       => [qw{ snps indels both any }],
            default     => q{both},
            store       => \$multiallelic_type,
            strict_type => 1,
        },
        outfile_path => {
            store       => \$outfile_path,
            strict_type => 1,
        },
        output_type => {
            allow       => [qw{ b u z v }],
            default     => q{v},
            store       => \$output_type,
            strict_type => 1,
        },
        reference_path => {
            defined     => 1,
            required    => 1,
            store       => \$reference_path,
            strict_type => 1,
        },
        regions_ref => { default => [], store => \$regions_ref, strict_type => 1, },
        samples_file_path => { store => \$samples_file_path, strict_type => 1, },
        samples_ref       => {
            default     => [],
            store       => \$samples_ref,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ bcftools norm };

    ## Bcftools base args
    @commands = bcftools_base(
        {
            commands_ref      => \@commands,
            regions_ref       => $regions_ref,
            outfile_path      => $outfile_path,
            output_type       => $output_type,
            samples_file_path => $samples_file_path,
            samples_ref       => $samples_ref,
        }
    );

    ## Options
    if ($multiallelic) {

        push @commands, q{--multiallelics} . $SPACE . $multiallelic . $multiallelic_type;
    }

    if ($reference_path) {

        push @commands, q{--fasta-ref} . $SPACE . $reference_path;
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
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );

    return @commands;
}

sub bcftools_reheader {

## Function : Perl wrapper for writing bcftools reheader recipe to already open $FILEHANDLE or return commands array. Based on bcftools 1.6.
## Returns  : @commands
## Arguments: $FILEHANDLE             => Filehandle to write to
##          : $infile_path            => Infile path
##          : $outfile_path           => Outfile path
##          : $regions_ref            => Regions to process {REF}
##          : $samples_file_path      => File of samples to annotate
##          : $samples_ref            => Samples to include or exclude if prefixed with "^"
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_path;
    my $regions_ref;
    my $samples_file_path;
    my $samples_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    my $tmpl = {
        FILEHANDLE   => { store   => \$FILEHANDLE, },
        infile_path  => { store   => \$infile_path, strict_type => 1, },
        outfile_path => { store   => \$outfile_path, strict_type => 1, },
        regions_ref  => { default => [], store => \$regions_ref, strict_type => 1, },
        samples_file_path => { store => \$samples_file_path, strict_type => 1, },
        samples_ref       => {
            default     => [],
            store       => \$samples_ref,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ bcftools reheader };

    ## Bcftools base args
    @commands = bcftools_base(
        {
            commands_ref => \@commands,
            regions_ref  => $regions_ref,
            outfile_path => $outfile_path,
        }
    );

    ## Options

    if ($samples_file_path) {

        push @commands, q{--samples} . $SPACE . $samples_file_path;
    }

    # Infile
    if ($infile_path) {

        push @commands, $infile_path;
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
            stdoutfile_path        => $stdoutfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );

    return @commands;
}

sub bcftools_rename_vcf_samples {

## Function : Rename vcf samples. The samples array will replace the sample names in the same order as supplied.
## Returns  :
## Arguments: $create_sample_file  => Create sample file for bcftools reheader
##          : $FILEHANDLE          => Filehandle to write to
##          : $index               => Generate index of reformated file
##          : $index_type          => Type of index
##          : $infile              => Vcf infile to rename samples for
##          : $outfile_path_prefix => Out file path no file_ending {Optional}
##          : $output_type         => Output type
##          : $sample_ids_ref      => Samples to rename in the same order as in the vcf {REF}
##          : $temp_directory      => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $create_sample_file;
    my $FILEHANDLE;
    my $infile;
    my $outfile_path_prefix;
    my $sample_ids_ref;
    my $temp_directory;

    ## Default(s)
    my $index;
    my $index_type;
    my $output_type;

    my $tmpl = {
        create_sample_file => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$create_sample_file,
            strict_type => 1,
        },
        FILEHANDLE => { defined => 1, required => 1, store => \$FILEHANDLE, },
        index      => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$index,
            strict_type => 1,
        },
        index_type => {
            allow       => [ undef, qw{ csi tbi } ],
            default     => q{csi},
            store       => \$index_type,
            strict_type => 1,
        },
        infile => { defined => 1, required => 1, store => \$infile, strict_type => 1, },
        outfile_path_prefix => { strict_type => 1, store => \$outfile_path_prefix, },
        output_type         => {
            allow       => [qw{ b u z v }],
            default     => q{v},
            store       => \$output_type,
            strict_type => 1,
        },
        sample_ids_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$sample_ids_ref,
            strict_type => 1,
        },
        temp_directory => {
            defined     => 1,
            required    => 1,
            store       => \$temp_directory,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    if ($create_sample_file) {

        bcftools_create_reheader_samples_file(
            {
                FILEHANDLE     => $FILEHANDLE,
                sample_ids_ref => $sample_ids_ref,
                temp_directory => $temp_directory,
            }
        );
    }

    ## Rename samples in VCF
    say {$FILEHANDLE} q{## Rename sample(s) names in VCF file};
    bcftools_reheader(
        {
            FILEHANDLE        => $FILEHANDLE,
            infile_path       => $infile,
            samples_file_path => catfile( $temp_directory, q{sample_name.txt} ),
        }
    );
    ## Pipe
    print {$FILEHANDLE} $PIPE . $SPACE;

    bcftools_view_and_index_vcf(
        {
            FILEHANDLE          => $FILEHANDLE,
            index               => $index,
            index_type          => $index_type,
            infile_path         => q{-},
            outfile_path_prefix => $outfile_path_prefix,
            output_type         => $output_type,
        }
    );

    return;
}

sub bcftools_create_reheader_samples_file {

## Function : Create reheader samples file.
## Returns  :
## Arguments: $FILEHANDLE          => Filehandle to write to
##          : $sample_ids_ref      => Samples to rename in the same order as in the vcf {REF}
##          : $temp_directory      => Temporary directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $sample_ids_ref;
    my $temp_directory;

    my $tmpl = {
        FILEHANDLE     => { defined => 1, required => 1, store => \$FILEHANDLE, },
        sample_ids_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$sample_ids_ref,
            strict_type => 1,
        },
        temp_directory => {
            defined     => 1,
            required    => 1,
            store       => \$temp_directory,
            strict_type => 1,
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
            FILEHANDLE      => $FILEHANDLE,
            format_string   => $format_string,
            stdoutfile_path => catfile( $temp_directory, q{sample_name.txt} ),
        }
    );
    say {$FILEHANDLE} $NEWLINE;
    return;
}

sub bcftools_roh {

## Function : Perl wrapper for writing bcftools roh recipe to $FILEHANDLE or return commands array. Based on bcftools 1.6.
## Returns  : @commands
## Arguments: $af_file_path           => Read allele frequencies from file (CHR\tPOS\tREF,ALT\tAF)
##          : $FILEHANDLE             => Filehandle to write to
##          : $infile_path            => Infile path to read from
##          : $outfile_path           => Outfile path to write to
##          : $output_type            => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##          : $regions_ref            => Regions to process {REF}
##          : $samples_file_path      => File of samples to annotate
##          : $samples_ref            => Samples to include or exclude if prefixed with "^"
##          : $skip_indels            => Skip indels as their genotypes are enriched for errors
##          : $stderrfile_path        => Stderr file path to write to
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile file path to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $af_file_path;
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_path;
    my $regions_ref;
    my $samples_file_path;
    my $samples_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;

    ## Default(s)
    my $skip_indels;
    my $output_type;

    my $tmpl = {
        af_file_path => { store => \$af_file_path, strict_type => 1, },
        FILEHANDLE   => { store => \$FILEHANDLE, },
        infile_path  => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        outfile_path => { store => \$outfile_path, strict_type => 1, },
        output_type  => {
            allow       => [qw{ b u z v}],
            default     => q{v},
            store       => \$output_type,
            strict_type => 1,
        },
        regions_ref => { default => [], store => \$regions_ref, strict_type => 1, },
        samples_file_path => { store => \$samples_file_path, strict_type => 1, },
        samples_ref       => {
            default     => [],
            store       => \$samples_ref,
            strict_type => 1,
        },
        samples_ref => {
            default     => [],
            store       => \$samples_ref,
            strict_type => 1,
        },
        skip_indels => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$skip_indels,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ bcftools roh };

    ## Bcftools base args
    @commands = bcftools_base(
        {
            commands_ref      => \@commands,
            regions_ref       => $regions_ref,
            output_type       => $output_type,
            outfile_path      => $outfile_path,
            samples_file_path => $samples_file_path,
            samples_ref       => $samples_ref,
        }
    );

    ## Options
    if ($af_file_path) {

        push @commands, q{--AF-file} . $SPACE . $af_file_path;
    }

    if ($skip_indels) {

        push @commands, q{--skip-indels};
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
            stdoutfile_path        => $stdoutfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );

    return @commands;
}

sub bcftools_stats {

## Function : Perl wrapper for writing bcftools stats recipe to already open $FILEHANDLE or return commands array. Based on bcftools 1.6.
## Returns  : @commands
## Arguments: $infile_path            => Infile path
##          : $outfile_path           => Outfile path
##          : $output_type            => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##          : $regions_ref            => Regions to process {REF}
##          : $samples_file_path      => File of samples to annotate
##          : $samples_ref            => Samples to include or exclude if prefixed with "^"
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $FILEHANDLE             => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $infile_path;
    my $outfile_path;
    my $regions_ref;
    my $samples_file_path;
    my $samples_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $FILEHANDLE;

    ## Default(s)
    my $output_type;

    my $tmpl = {
        FILEHANDLE   => { store => \$FILEHANDLE, },
        infile_path  => { store => \$infile_path, strict_type => 1, },
        outfile_path => { store => \$outfile_path, strict_type => 1, },
        output_type => {
            allow       => [qw{ b u z v}],
            default     => q{v},
            store       => \$output_type,
            strict_type => 1,
        },
        regions_ref => { default => [], store => \$regions_ref, strict_type => 1, },
        samples_file_path => { store => \$samples_file_path, strict_type => 1, },
        samples_ref       => {
            default     => [],
            store       => \$samples_ref,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ bcftools stats };

    ## Bcftools base args
    @commands = bcftools_base(
        {
            commands_ref      => \@commands,
            regions_ref       => $regions_ref,
            samples_file_path => $samples_file_path,
            samples_ref       => $samples_ref,
        }
    );

    ## Infile
    if ($infile_path) {

        push @commands, $infile_path;
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
            stdoutfile_path        => $stdoutfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub bcftools_view {

## Function : Perl wrapper for writing bcftools view recipe to $FILEHANDLE or return commands array. Based on bcftools 1.6.
## Returns  : @commands
## Arguments: $apply_filters_ref      => Require at least one of the listed FILTER strings
##          : $exclude_types_ref      => Exclude comma-separated list of variant types: snps,indels,mnps,other
##          : $exclude                => Exclude sites for which the expression is true
##          : $FILEHANDLE             => Filehandle to write to
##          : $genotype               => Genotype to include (hom|het|miss). Prefix with "^" for exclude
##          : $include                => Include only sites for which the expression is true
##          : $infile_path            => Infile path to read from
##          : $max_alleles            => Max alleles listed in REF and ALT columns
##          : $min_alleles            => Min alleles listed in REF and ALT columns
##          : $outfile_path           => Outfile path to write to
##          : $output_type            => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
##          : $regions_ref            => Regions to process {REF}
##          : $samples_file_path      => File of samples to annotate
##          : $samples_ref            => Samples to include or exclude if prefixed with "^"
##          : $stderrfile_path        => Stderr file path to write to
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile file path to write to
##          : $types                  => Comma separated variant types to include (snps|indels|mnps|other), based on based on REF,ALT

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $apply_filters_ref;
    my $exclude_types_ref;
    my $exclude;
    my $FILEHANDLE;
    my $genotype;
    my $include;
    my $infile_path;
    my $max_alleles;
    my $min_alleles;
    my $outfile_path;
    my $regions_ref;
    my $samples_file_path;
    my $samples_ref;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $types;

    ## Default(s)
    my $output_type;

    my $tmpl = {
        FILEHANDLE => { store => \$FILEHANDLE, },
        apply_filters_ref =>
          { default => [], store => \$apply_filters_ref, strict_type => 1, },
        exclude_types_ref =>
          { default => [], store => \$exclude_types_ref, strict_type => 1, },
        exclude  => { store => \$exclude, strict_type => 1, },
        genotype => {
            store       => \$genotype,
            strict_type => 1,
        },
        include     => { store => \$include,     strict_type => 1, },
        infile_path => { store => \$infile_path, strict_type => 1, },
        max_alleles => {
            allow       => [ undef, qr/ ^\d+$ /xms ],
            store       => \$max_alleles,
            strict_type => 1,
        },
        min_alleles => {
            allow       => [ undef, qr/ ^\d+$ /xms ],
            store       => \$min_alleles,
            strict_type => 1,
        },
        outfile_path => { store => \$outfile_path, strict_type => 1, },
        output_type  => {
            allow       => [qw{ b u z v}],
            default     => q{v},
            store       => \$output_type,
            strict_type => 1,
        },
        regions_ref => { default => [], store => \$regions_ref, strict_type => 1, },
        samples_file_path => { store => \$samples_file_path, strict_type => 1, },
        samples_ref       => {
            default     => [],
            store       => \$samples_ref,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        stderrfile_path_append =>
          { store => \$stderrfile_path_append, strict_type => 1, },
        stdoutfile_path => { store => \$stdoutfile_path, strict_type => 1, },
        types           => {
            store       => \$types,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ bcftools view };

    ## Bcftools base args
    @commands = bcftools_base(
        {
            commands_ref      => \@commands,
            regions_ref       => $regions_ref,
            output_type       => $output_type,
            samples_file_path => $samples_file_path,
            samples_ref       => $samples_ref,
        }
    );

    ## Options
    if ( @{$apply_filters_ref} ) {

        push @commands, q{--apply-filters} . $SPACE . join $COMMA, @{$apply_filters_ref};
    }

    if ( @{$exclude_types_ref} ) {

        push @commands, q{--exclude-types} . $SPACE . join $COMMA, @{$exclude_types_ref};
    }

    if ($exclude) {

        push @commands, q{--exclude} . $SPACE . $exclude;
    }

    if ($genotype) {

        push @commands, q{--genotype} . $SPACE . $genotype;
    }

    if ($include) {

        push @commands, q{--include} . $SPACE . $include;
    }

    if ($max_alleles) {

        push @commands, q{--max-alleles} . $SPACE . $max_alleles;
    }

    if ($min_alleles) {

        push @commands, q{--min-alleles} . $SPACE . $min_alleles;
    }

    if ($outfile_path) {

        push @commands, q{--output-file} . $SPACE . $outfile_path;
    }

    if ($types) {

        push @commands, q{--types} . $SPACE . $types;
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
            stdoutfile_path        => $stdoutfile_path,
        }
      );

    unix_write_to_file(
        {
            commands_ref => \@commands,
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,
        }
    );

    return @commands;
}

sub bcftools_view_and_index_vcf {

## Function : View variant calling file and index.
## Returns  :
## Arguments: $FILEHANDLE          => SBATCH script FILEHANDLE to print to
##          : $index               => Generate index of reformated file
##          : $index_type          => Type of index
##          : $infile_path         => Path to infile to compress and index
##          : $outfile_path_prefix => Out file path no file_ending {Optional}
##          : $output_type         => 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $infile_path;
    my $outfile_path_prefix;

    ## Default(s)
    my $index;
    my $index_type;
    my $output_type;

    my $tmpl = {
        FILEHANDLE  => { defined => 1, required => 1, store => \$FILEHANDLE, },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
        index => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$index,
            strict_type => 1,
        },
        index_type => {
            allow       => [ undef, qw{ csi tbi } ],
            default     => q{csi},
            store       => \$index_type,
            strict_type => 1,
        },
        outfile_path_prefix => { store => \$outfile_path_prefix, strict_type => 1, },
        output_type         => {
            allow       => [qw{ b u z v }],
            default     => q{v},
            store       => \$output_type,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $outfile_path;
    my %output_type_ending = (
        b => q{.bcf},
        u => q{.bcf},
        v => q{.vcf},
        z => q{.vcf.gz},
    );

    if ( defined $outfile_path_prefix ) {

        $outfile_path = $outfile_path_prefix . $output_type_ending{$output_type};
    }

    bcftools_view(
        {
            FILEHANDLE   => $FILEHANDLE,
            infile_path  => $infile_path,
            outfile_path => $outfile_path,
            output_type  => $output_type,
        }
    );
    say {$FILEHANDLE} $NEWLINE;

    if ($index) {

        say {$FILEHANDLE} q{## Index};

        bcftools_index(
            {
                FILEHANDLE  => $FILEHANDLE,
                infile_path => $outfile_path,
                output_type => $index_type,
            }
        );
        say {$FILEHANDLE} $NEWLINE;
    }
    return;
}

1;
