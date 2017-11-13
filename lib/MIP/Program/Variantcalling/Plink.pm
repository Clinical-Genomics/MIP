package MIP::Program::Variantcalling::Plink;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;

## CPANM
use Readonly;

## MIPs lib/
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ plink plink_variant_pruning plink_fix_fam_ped_map_freq plink_calculate_inbreeding plink_create_mibs plink_check_sex_chroms plink_sex_check };
}

## Constants
Readonly my $SPACE => q{ };
Readonly my $COMMA => q{,};

sub plink_variant_pruning {

## Function : Perl wrapper for writing Plink recipe to prune variants and create unique IDs. Based on Plink 1.90p 64-bit (25 Mar 2016).
## Returns  : @commands
## Arguments: $vcffile_path        => Vcf file path
##          : $outfile_prefix      => Outfile prefix
##          : $stderrfile_path     => Stderrfile path
##          : $vcf_require_gt      => Skip variants where the GT field is absent
##          : $vcf_half_call       => Specify how '0/.' and similar VCF GT values should be handled
##          : $set_missing_var_ids => Assign chromosome-and-position-based IDs
##          : $const_fid           => Assign family id to all samples in file
##          : $make_bed            => Save data in plink binary format
##          : $indep               => Produce a pruned subset of markers that are in approximate linkage equilibrium with each other
##          : $indep_window_size   => Indep window size (kb)
##          : $indep_step_size     => Indep step size (variant ct)
##          : $indep_vif_threshold => Indep vif threshold
##          : $FILEHANDLE          => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $outfile_prefix;
    my $stderrfile_path;
    my $vcffile_path;
    my $vcf_half_call;
    my $set_missing_var_ids;
    my $const_fid;
    my $make_bed;
    my $indep;
    my $indep_window_size;
    my $indep_step_size;
    my $indep_vif_threshold;
    my $FILEHANDLE;

    ## Default(s)
    my $vcf_require_gt;

    my $tmpl = {
        vcffile_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$vcffile_path,
        },
        outfile_prefix => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$outfile_prefix,
        },
        vcf_require_gt => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$vcf_require_gt,
        },
        vcf_half_call => {
            allow       => [ undef, qw{ error haploid missing reference } ],
            strict_type => 1,
            store => \$vcf_half_call,
        },
        set_missing_var_ids => {
            required    => 1,
            strict_type => 1,
            store       => \$set_missing_var_ids,
        },
        const_fid => {
            required    => 1,
            strict_type => 1,
            store       => \$const_fid
        },
        make_bed => {
            strict_type => 1,
            store       => \$make_bed,
        },
        indep => {
            required    => 1,
            strict_type => 1,
            store       => \$indep,
        },
        indep_window_size =>
          { required => 1, strict_type => 1, store => \$indep_window_size, },
        indep_step_size =>
          { required => 1, strict_type => 1, store => \$indep_step_size, },
        indep_vif_threshold =>
          { required => 1, strict_type => 1, store => \$indep_vif_threshold, },
        FILEHANDLE      => { store       => \$FILEHANDLE, },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },

    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = q{plink2};

    push @commands, q{--vcf} . $SPACE . $vcffile_path;

    if ($vcf_require_gt) {

        push @commands, q{--vcf-require-gt};
    }

    if ($vcf_half_call) {

        push @commands, q{--vcf-half-call} . $SPACE . $vcf_half_call;
    }

    push @commands, q{--set-missing-var-ids} . $SPACE . $set_missing_var_ids;

    push @commands, q{--const-fid} . $SPACE . $const_fid;

    if ($make_bed) {

        push @commands, q{--make-bed};
    }

    push @commands, q{--indep};
    push @commands, $indep_window_size;
    push @commands, $indep_step_size;
    push @commands, $indep_vif_threshold;

    ## Outfile
    push @commands, q{--out} . $SPACE . $outfile_prefix;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
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

sub plink_fix_fam_ped_map_freq {

## Function : Perl wrapper for writing Plink recipe to create ped, map, frequency report. Updates fam file as well. Based on Plink 1.90p 64-bit (25 Mar 2016).
## Returns  : @commands
## Arguments: $binary_fileset_prefix => Specify .bed + .bim + .fam prefix
##          : $fam_file_path         => Fam file path
##          : $make_just_fam         => Just update fam file
##          : $recode                => Create a new text fileset with all filters applied
##          : $freqx                 => Writes a minor allele frequency report
##          : $outfile_prefix        => Outfile prefix
##          : $FILEHANDLE            => Filehandle to write to
##          : $stderrfile_path       => Stderrfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $binary_fileset_prefix;
    my $fam_file_path;
    my $make_just_fam;
    my $outfile_prefix;
    my $FILEHANDLE;
    my $stderrfile_path;

    ## Default(s)
    my $recode;
    my $freqx;

    my $tmpl = {
        binary_fileset_prefix => {
            required    => 1,
            strict_type => 1,
            store       => \$binary_fileset_prefix,
        },
        fam_file_path => {
            required    => 1,
            strict_type => 1,
            store       => \$fam_file_path,
        },
        make_just_fam => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$make_just_fam,
        },
        recode => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$recode,
        },
        freqx => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$freqx,
        },
        outfile_prefix => {
            required    => 1,
            strict_type => 1,
            store       => \$outfile_prefix,
        },
        FILEHANDLE      => { store       => \$FILEHANDLE },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },

    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = q{plink2};

    push @commands, q{--bfile} . $SPACE . $binary_fileset_prefix;

    push @commands, q{--fam} . $SPACE . $fam_file_path;

    if ($make_just_fam) {

        push @commands, q{--make-just-fam};
    }

    if ($recode) {

        push @commands, q{--recode};
    }

    if ($freqx) {

        push @commands, q{--freqx};
    }

    push @commands, q{--out} . $SPACE . $outfile_prefix;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
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

sub plink_calculate_inbreeding {

## Function : Perl wrapper for writing Plink recipe to calculate inbreeding coefficients per family. Based on Plink 1.90p 64-bit (25 Mar 2016).
## Returns  : @commands
## Arguments: $binary_fileset_prefix   => Specify .bed + .bim + .fam prefix
##          : $outfile_prefix          => Outfile prefix
##          : $het                     => Reports method-of-moments estimates
##          : $small_sample            => Inlcude n/(n-1) multiplier in Nei's expected homozygosity formula
##          : $inbreeding_coefficients => Estimate inbreeding coefficients
##          : $extract_file            => Exclude all variants not named in the file
##          : $FILEHANDLE              => Filehandle to write to
##          : $stderrfile_path         => Stderrfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $binary_fileset_prefix;
    my $outfile_prefix;
    my $stderrfile_path;
    my $inbreeding_coefficients;
    my $extract_file;
    my $FILEHANDLE;

    ## Default(s)
    my $het;
    my $small_sample;

    my $tmpl = {
        binary_fileset_prefix => {
            required    => 1,
            strict_type => 1,
            store       => \$binary_fileset_prefix
        },
        outfile_prefix => {
            required    => 1,
            strict_type => 1,
            store       => \$outfile_prefix,
        },
        het => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$het
        },
        small_sample => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$small_sample
        },
        inbreeding_coefficients => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$inbreeding_coefficients
        },
        extract_file => {
            required    => 1,
            strict_type => 1,
            store       => \$extract_file
        },
        FILEHANDLE      => { store       => \$FILEHANDLE },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },

    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = q{plink2};

    push @commands, q{--bfile} . $SPACE . $binary_fileset_prefix;

    if ($het) {

        push @commands, q{--het};
    }

    if ($small_sample) {

        push @commands, q{small-sample};
    }

    if ($inbreeding_coefficients) {

        push @commands, q{--ibc};
    }

    push @commands, q{--extract} . $SPACE . $extract_file;

    push @commands, q{--out} . $SPACE . $outfile_prefix;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
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

sub plink_create_mibs {

## Function : Perl wrapper for writing Plink recipe to create .mibs per family. Based on Plink 1.90p 64-bit (25 Mar 2016).
## Returns  : @commands
## Arguments: $ped_file_path   => Ped file path
##          : $outfile_prefix  => Outfile prefix
##          : $map_file_path   => Map file path
##          : $cluster         => Perform IBS clustering
##          : $matrix          => Create a N x N matrix of genome-wide average IBS pairwise identities
##          : $FILEHANDLE      => Filehandle to write to
##          : $stderrfile_path => Stderrfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $ped_file_path;
    my $outfile_prefix;
    my $map_file_path;
    my $FILEHANDLE;
    my $stderrfile_path;

    ## Default(s)
    my $cluster;
    my $matrix;

    my $tmpl = {
        ped_file_path =>
          { required => 1, strict_type => 1, store => \$ped_file_path },
        map_file_path =>
          { required => 1, strict_type => 1, store => \$map_file_path },
        cluster => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$cluster
        },
        matrix => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$matrix
        },
        outfile_prefix =>
          { required => 1, strict_type => 1, store => \$outfile_prefix },
        FILEHANDLE      => { store       => \$FILEHANDLE },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = q{plink2};

    push @commands, q{--ped} . $SPACE . $ped_file_path;

    push @commands, q{--map} . $SPACE . $map_file_path;

    if ($cluster) {

        push @commands, q{--cluster};
    }

    if ($matrix) {

        push @commands, q{--matrix};
    }

    push @commands, q{--out} . $SPACE . $outfile_prefix;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
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

sub plink_check_sex_chroms {

## Function : Perl wrapper for writing Plink recipe to check sex chromosomes. Based on Plink 1.90p 64-bit (25 Mar 2016).
## Returns  : @commands
## Arguments: $regions_ref   => The regions to process {REF}
##          : $split_x       => Changes the chromosome codes of all variants in the region to XY
##          : $no_fail                 => Do not fail when no variants would be affected by split_x
##          : $make_bed                => Save data in plink binary format
##          : $binary_fileset_prefix   => Specify .bed + .bim + .fam prefix
##          : $outfile_prefix          => Outfile prefix
##          : $FILEHANDLE              => Filehandle to write to
##          : $stderrfile_path         => Stderrfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $regions_ref;
    my $split_x;
    my $make_bed;
    my $binary_fileset_prefix;
    my $outfile_prefix;
    my $FILEHANDLE;
    my $stderrfile_path;

    ## Default(s)
    my $no_fail;

    my $tmpl = {
        regions_ref => {
            required    => 1,
            default     => [],
            strict_type => 1,
            store       => \$regions_ref
        },
        split_x => { required => 1, strict_type => 1, store => \$split_x },
        no_fail => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$no_fail,
        },
        make_bed => { strict_type => 1, store => \$make_bed },
        binary_fileset_prefix =>
          { required => 1, strict_type => 1, store => \$binary_fileset_prefix },
        outfile_prefix =>
          { required => 1, strict_type => 1, store => \$outfile_prefix },
        FILEHANDLE      => { store       => \$FILEHANDLE },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = q{plink2};

    push @commands, q{--chr} . $SPACE . join $COMMA, @{$regions_ref};

    push @commands, q{--split-x} . $SPACE . $split_x;

    if ($no_fail) {

        push @commands, q{no-fail};
    }

    if ($make_bed) {

        push @commands, q{--make-bed};
    }

    push @commands, q{--bfile} . $SPACE . $binary_fileset_prefix;

    push @commands, q{--out} . $SPACE . $outfile_prefix;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
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

sub plink_sex_check {

## Function : Perl wrapper for writing Plink recipe to do sex check. Based on Plink 1.90p 64-bit (25 Mar 2016).
## Returns  : @commands
## Arguments: $sex_check_min_f   => Sex check minimum F [female male]
##          : $extract_file            => Exclude all variants not named in the file
##          : $read_freqfile_path      => Read from frequency file path
##          : $binary_fileset_prefix   => Specify .bed + .bim + .fam prefix
##          : $outfile_prefix          => Outfile prefix
##          : $FILEHANDLE              => Filehandle to write to
##          : $stderrfile_path         => Stderrfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $sex_check_min_f;
    my $extract_file;
    my $read_freqfile_path;
    my $binary_fileset_prefix;
    my $outfile_prefix;
    my $stderrfile_path;
    my $FILEHANDLE;

    my $tmpl = {
        sex_check_min_f =>
          { required => 1, strict_type => 1, store => \$sex_check_min_f },
        extract_file =>
          { required => 1, strict_type => 1, store => \$extract_file },
        read_freqfile_path =>
          { required => 1, strict_type => 1, store => \$read_freqfile_path },
        binary_fileset_prefix =>
          { required => 1, strict_type => 1, store => \$binary_fileset_prefix },
        outfile_prefix =>
          { required => 1, strict_type => 1, store => \$outfile_prefix },
        FILEHANDLE      => { store       => \$FILEHANDLE },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },

    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = q{plink2 --check-sex};

    push @commands, $sex_check_min_f;

    push @commands, q{--extract} . $SPACE . $extract_file;

    push @commands, q{--read-freq} . $SPACE . $read_freqfile_path;

    push @commands, q{--bfile} . $SPACE . $binary_fileset_prefix;

    push @commands, q{--out} . $SPACE . $outfile_prefix;

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path => $stderrfile_path,
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

1;
