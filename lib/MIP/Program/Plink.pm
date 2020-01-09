package MIP::Program::Plink;

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
use MIP::Constants qw{ $COMMA $SPACE };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ plink_calculate_inbreeding
      plink_check_sex_chroms
      plink_create_mibs
      plink_fix_fam_ped_map_freq
      plink_sex_check plink_variant_pruning };
}

sub plink_variant_pruning {

## Function : Perl wrapper for writing Plink recipe to prune variants and create unique IDs. Based on Plink 1.90p 64-bit (25 Mar 2016).
## Returns  : @commands
## Arguments: $const_fid           => Assign case id to all samples in file
##          : $filehandle          => Filehandle to write to
##          : $indep               => Produce a pruned subset of markers that are in approximate linkage equilibrium with each other
##          : $indep_step_size     => Indep step size (variant ct)
##          : $indep_vif_threshold => Indep vif threshold
##          : $indep_window_size   => Indep window size (kb)
##          : $make_bed            => Save data in plink binary format
##          : $outfile_prefix      => Outfile prefix
##          : $set_missing_var_ids => Assign chromosome-and-position-based IDs
##          : $stderrfile_path     => Stderrfile path
##          : $vcf_half_call       => Specify how '0/.' and similar VCF GT values should be handled
##          : $vcf_require_gt      => Skip variants where the GT field is absent
##          : $vcffile_path        => Vcf file path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $const_fid;
    my $filehandle;
    my $indep;
    my $indep_step_size;
    my $indep_vif_threshold;
    my $indep_window_size;
    my $make_bed;
    my $outfile_prefix;
    my $set_missing_var_ids;
    my $stderrfile_path;
    my $vcffile_path;
    my $vcf_half_call;

    ## Default(s)
    my $vcf_require_gt;

    my $tmpl = {
        const_fid => {
            required    => 1,
            store       => \$const_fid,
            strict_type => 1,
        },
        filehandle => { store => \$filehandle, },
        indep      => {
            required    => 1,
            store       => \$indep,
            strict_type => 1,
        },
        indep_step_size =>
          { required => 1, store => \$indep_step_size, strict_type => 1, },
        indep_vif_threshold =>
          { required => 1, store => \$indep_vif_threshold, strict_type => 1, },
        indep_window_size =>
          { required => 1, store => \$indep_window_size, strict_type => 1, },
        make_bed => {
            store       => \$make_bed,
            strict_type => 1,
        },
        outfile_prefix => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_prefix,
            strict_type => 1,
        },
        set_missing_var_ids => {
            required    => 1,
            store       => \$set_missing_var_ids,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
        vcffile_path    => {
            defined     => 1,
            required    => 1,
            store       => \$vcffile_path,
            strict_type => 1,
        },
        vcf_half_call => {
            allow       => [ undef, qw{ error haploid missing reference } ],
            store       => \$vcf_half_call,
            strict_type => 1,
        },
        vcf_require_gt => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$vcf_require_gt,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ plink2 };

    if ($make_bed) {

        push @commands, q{--make-bed};
    }

    if ($vcf_require_gt) {

        push @commands, q{--vcf-require-gt};
    }

    if ($vcf_half_call) {

        push @commands, q{--vcf-half-call} . $SPACE . $vcf_half_call;
    }

    push @commands, q{--set-missing-var-ids} . $SPACE . $set_missing_var_ids;

    push @commands, q{--indep};
    push @commands, $indep_window_size;
    push @commands, $indep_step_size;
    push @commands, $indep_vif_threshold;

    push @commands, q{--vcf} . $SPACE . $vcffile_path;

    push @commands, q{--const-fid} . $SPACE . $const_fid;

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
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub plink_fix_fam_ped_map_freq {

## Function : Perl wrapper for writing Plink recipe to create ped, map, frequency report. Updates fam file as well. Based on Plink 1.90p 64-bit (25 Mar 2016).
## Returns  : @commands
## Arguments: $allow_no_sex          => Allow no sex for sample ids
##          : $binary_fileset_prefix => Specify .bed + .bim + .fam prefix
##          : $fam_file_path         => Fam file path
##          : $filehandle            => Filehandle to write to
##          : $freqx                 => Writes a minor allele frequency report
##          : $make_just_fam         => Just update fam file
##          : $outfile_prefix        => Outfile prefix
##          : $recode                => Create a new text fileset with all filters applied
##          : $stderrfile_path       => Stderrfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $allow_no_sex;
    my $binary_fileset_prefix;
    my $fam_file_path;
    my $filehandle;
    my $make_just_fam;
    my $outfile_prefix;
    my $stderrfile_path;

    ## Default(s)
    my $freqx;
    my $recode;

    my $tmpl = {
        allow_no_sex => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$allow_no_sex,
            strict_type => 1,
        },
        binary_fileset_prefix => {
            required    => 1,
            store       => \$binary_fileset_prefix,
            strict_type => 1,
        },
        fam_file_path => {
            required    => 1,
            store       => \$fam_file_path,
            strict_type => 1,
        },
        filehandle => { store => \$filehandle, },
        freqx      => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$freqx,
            strict_type => 1,
        },
        make_just_fam => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$make_just_fam,
            strict_type => 1,
        },
        outfile_prefix => {
            required    => 1,
            store       => \$outfile_prefix,
            strict_type => 1,
        },
        recode => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$recode,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ plink2 };

    push @commands, q{--bfile} . $SPACE . $binary_fileset_prefix;

    if ($make_just_fam) {

        push @commands, q{--make-just-fam};
    }

    push @commands, q{--fam} . $SPACE . $fam_file_path;

    if ($recode) {

        push @commands, q{--recode};
    }

    if ($freqx) {

        push @commands, q{--freqx};
    }
    if ($allow_no_sex) {

        push @commands, q{--allow-no-sex};
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
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub plink_calculate_inbreeding {

## Function : Perl wrapper for writing Plink recipe to calculate inbreeding coefficients per case. Based on Plink 1.90p 64-bit (25 Mar 2016).
## Returns  : @commands
## Arguments: $binary_fileset_prefix   => Specify .bed + .bim + .fam prefix
##          : $extract_file            => Exclude all variants not named in the file
##          : $filehandle              => Filehandle to write to
##          : $het                     => Reports method-of-moments estimates
##          : $inbreeding_coefficients => Estimate inbreeding coefficients
##          : $outfile_prefix          => Outfile prefix
##          : $small_sample            => Inlcude n/(n-1) multiplier in Nei's expected homozygosity formula
##          : $stderrfile_path         => Stderrfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $binary_fileset_prefix;
    my $extract_file;
    my $filehandle;
    my $inbreeding_coefficients;
    my $outfile_prefix;
    my $stderrfile_path;

    ## Default(s)
    my $het;
    my $small_sample;

    my $tmpl = {
        binary_fileset_prefix => {
            required    => 1,
            store       => \$binary_fileset_prefix,
            strict_type => 1,
        },
        extract_file => {
            required    => 1,
            store       => \$extract_file,
            strict_type => 1,
        },
        filehandle => { store => \$filehandle, },
        het        => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$het,
            strict_type => 1,
        },
        inbreeding_coefficients => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$inbreeding_coefficients,
            strict_type => 1,
        },
        outfile_prefix => {
            required    => 1,
            store       => \$outfile_prefix,
            strict_type => 1,
        },
        small_sample => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$small_sample,
            strict_type => 1,
        },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },

    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ plink2 };

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
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub plink_create_mibs {

## Function : Perl wrapper for writing Plink recipe to create .mibs per case. Based on Plink 1.90p 64-bit (25 Mar 2016).
## Returns  : @commands
## Arguments: $cluster         => Perform IBS clustering
##          : $filehandle      => Filehandle to write to
##          : $map_file_path   => Map file path
##          : $matrix          => Create a N x N matrix of genome-wide average IBS pairwise identities
##          : $outfile_prefix  => Outfile prefix
##          : $ped_file_path   => Ped file path
##          : $stderrfile_path => Stderrfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $map_file_path;
    my $outfile_prefix;
    my $ped_file_path;
    my $stderrfile_path;

    ## Default(s)
    my $cluster;
    my $matrix;

    my $tmpl = {
        cluster => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$cluster,
            strict_type => 1,
        },
        filehandle    => { store    => \$filehandle, },
        map_file_path => { required => 1, store => \$map_file_path, strict_type => 1, },
        matrix => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$matrix,
            strict_type => 1,
        },
        outfile_prefix => { required => 1, store => \$outfile_prefix, strict_type => 1, },
        ped_file_path  => { required => 1, store => \$ped_file_path,  strict_type => 1, },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ plink2 };

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
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub plink_check_sex_chroms {

## Function : Perl wrapper for writing Plink recipe to check sex chromosomes. Based on Plink 1.90p 64-bit (25 Mar 2016).
## Returns  : @commands
## Arguments: $binary_fileset_prefix => Specify .bed + .bim + .fam prefix
##          : $filehandle            => Filehandle to write to
##          : $make_bed              => Save data in plink binary format
##          : $no_fail               => Do not fail when no variants would be affected by split_x
##          : $outfile_prefix        => Outfile prefix
##          : $regions_ref           => The regions to process {REF}
##          : $split_x               => Changes the chromosome codes of all variants in the region to XY
##          : $stderrfile_path       => Stderrfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $binary_fileset_prefix;
    my $filehandle;
    my $make_bed;
    my $outfile_prefix;
    my $regions_ref;
    my $split_x;
    my $stderrfile_path;

    ## Default(s)
    my $no_fail;

    my $tmpl = {
        binary_fileset_prefix => {
            required    => 1,
            store       => \$binary_fileset_prefix,
            strict_type => 1,
        },
        filehandle => { store => \$filehandle, },
        make_bed   => { store => \$make_bed, strict_type => 1, },
        no_fail => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$no_fail,
            strict_type => 1,
        },
        outfile_prefix => { required => 1, store => \$outfile_prefix, strict_type => 1, },
        regions_ref    => {
            default     => [],
            required    => 1,
            store       => \$regions_ref,
            strict_type => 1,
        },
        split_x => { required => 1, store => \$split_x, strict_type => 1, },
        stderrfile_path => { store => \$stderrfile_path, strict_type => 1, },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ plink2 };

    push @commands, q{--bfile} . $SPACE . $binary_fileset_prefix;

    push @commands, q{--split-x} . $SPACE . $split_x;

    if ($no_fail) {

        push @commands, q{no-fail};
    }

    if ($make_bed) {

        push @commands, q{--make-bed};
    }

    push @commands, q{--chr} . $SPACE . join $COMMA, @{$regions_ref};

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
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

sub plink_sex_check {

## Function : Perl wrapper for writing Plink recipe to do sex check. Based on Plink 1.90p 64-bit (25 Mar 2016).
## Returns  : @commands
## Arguments: $binary_fileset_prefix => Specify .bed + .bim + .fam prefix
##          : $extract_file          => Exclude all variants not named in the file
##          : $filehandle            => Filehandle to write to
##          : $outfile_prefix        => Outfile prefix
##          : $read_freqfile_path    => Read from frequency file path
##          : $sex_check_min_f       => Sex check minimum F [female male]
##          : $stderrfile_path       => Stderrfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $binary_fileset_prefix;
    my $extract_file;
    my $filehandle;
    my $outfile_prefix;
    my $read_freqfile_path;
    my $sex_check_min_f;
    my $stderrfile_path;

    my $tmpl = {
        binary_fileset_prefix => {
            required    => 1,
            store       => \$binary_fileset_prefix,
            strict_type => 1,
        },
        extract_file   => { store    => \$extract_file, strict_type => 1, },
        filehandle     => { store    => \$filehandle, },
        outfile_prefix => { required => 1, store => \$outfile_prefix, strict_type => 1, },
        read_freqfile_path => { store => \$read_freqfile_path, strict_type => 1, },
        sex_check_min_f    => { store => \$sex_check_min_f, strict_type => 1, },
        stderrfile_path    => { store => \$stderrfile_path, strict_type => 1, },

    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ plink2 };

    push @commands, q{--bfile} . $SPACE . $binary_fileset_prefix;

    push @commands, q{--check-sex};

    if ($sex_check_min_f) {

        push @commands, $sex_check_min_f;
    }

    if ($read_freqfile_path) {

        push @commands, q{--read-freq} . $SPACE . $read_freqfile_path;
    }

    if ($extract_file) {

        push @commands, q{--extract} . $SPACE . $extract_file;
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
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );
    return @commands;
}

1;
