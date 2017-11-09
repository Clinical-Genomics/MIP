package MIP::Program::Variantcalling::Plink;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };

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
    our @EXPORT_OK = qw{ plink };
}

## Constants
Readonly my $SPACE => q{ };
Readonly my $COMMA => q{,};

sub plink {

## Function : Perl wrapper for writing Plink recipe to already open $FILEHANDLE or return commands array. Based on Plink 1.90p 64-bit (25 Mar 2016).
## Returns  : @commands
## Arguments:
##          : $regions_ref             => The regions to process {REF}
##          : $vcffile_path            => Vcf file path
##          : $outfile_prefix          => Outfile prefix
##          : $stderrfile_path         => Stderrfile path
##          : $stdoutfile_path         => Stdoutfile path
##          : $FILEHANDLE              => Filehandle to write to
##          : $make_bed                => Save data in plink binary format
##          : $binary_fileset_prefix   => Specify .bed + .bim + .fam prefix
##          : $split_x                 => Changes the chromosome codes of all variants in the region to XY
##          : $set_missing_var_ids     => Assign chromosome-and-position-based IDs
##          : $indep                   => Produce a pruned subset of markers that are in approximate linkage equilibrium with each other
##          : $indep_window_size       => Indep window size (kb)
##          : $indep_step_size         => Indep step size (variant ct)
##          : $indep_vif_threshold     => Indep vif threshold
##          : $extract_file            => Exclude all variants not named in the file
##          : $fam_file_path           => Fam file path
##          : $ped_file_path           => Ped file path
##          : $map_file_path           => Map file path
##          : $keep_file_path          => Keep file path
##          : $read_freqfile_path      => Read from frequency file path
##          : $const_fid               => Assign family id to all samples in file
##          : $sex_check_min_f         => Sex check minimum F [female male]
##          : $vcf_half_call           => Specify how '0/.' and similar VCF GT values should be handled
##          : $no_fail                 => Do not fail when no variants would be affected by split_x
##          : $check_sex               => Normally compares sex assignments in the input dataset with those imputed from X chromosome inbreeding coefficients
##          : $y_only                  => {female max Y obs} {male min Y obs}
##          : $vcf_require_gt          => Skip variants where the GT field is absent
##          : $het                     => Reports method-of-moments estimates
##          : $small_sample            => Inlcude n/(n-1) multiplier in Nei's expected homozygosity formula
##          : $inbreeding_coefficients => Estimate inbreeding coefficients
##          : $recode                  => Create a new text fileset with all filters applied
##          : $make_just_fam           => Just update fam file
##          : $cluster                 => Perform IBS clustering
##          : $matrix                  => Create a N x N matrix of genome-wide average IBS pairwise identities
##          : $freqx                   => Writes a minor allele frequency report
##          : $allow_no_sex            => Allow no sex for sample ids

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $regions_ref;
    my $vcffile_path;
    my $outfile_prefix;
    my $stderrfile_path;
    my $stdoutfile_path;
    my $FILEHANDLE;
    my $make_bed;
    my $binary_fileset_prefix;
    my $split_x;
    my $set_missing_var_ids;
    my $indep;
    my $indep_window_size;
    my $indep_step_size;
    my $indep_vif_threshold;
    my $extract_file;
    my $fam_file_path;
    my $ped_file_path;
    my $map_file_path;
    my $keep_file_path;
    my $read_freqfile_path;
    my $const_fid;
    my $sex_check_min_f;
    my $vcf_half_call;
    my $inbreeding_coefficients;
    my $allow_no_sex;

    ## Default(s)
    my $no_fail;
    my $check_sex;
    my $y_only;
    my $vcf_require_gt;
    my $het;
    my $small_sample;
    my $recode;
    my $make_just_fam;
    my $cluster;
    my $matrix;
    my $freqx;

    my $tmpl = {
        regions_ref =>
          { default => [], strict_type => 1, store => \$regions_ref },
        vcffile_path    => { strict_type => 1, store => \$vcffile_path },
        outfile_prefix  => { strict_type => 1, store => \$outfile_prefix },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
        FILEHANDLE => { store       => \$FILEHANDLE },
        make_bed   => { strict_type => 1, store => \$make_bed },
        binary_fileset_prefix =>
          { strict_type => 1, store => \$binary_fileset_prefix },
        split_x => { strict_type => 1, store => \$split_x },
        set_missing_var_ids =>
          { strict_type => 1, store => \$set_missing_var_ids },
        indep             => { strict_type => 1, store => \$indep },
        indep_window_size => { strict_type => 1, store => \$indep_window_size },
        indep_step_size   => { strict_type => 1, store => \$indep_step_size },
        indep_vif_threshold =>
          { strict_type => 1, store => \$indep_vif_threshold },
        extract_file   => { strict_type => 1, store => \$extract_file },
        fam_file_path  => { strict_type => 1, store => \$fam_file_path },
        ped_file_path  => { strict_type => 1, store => \$ped_file_path },
        map_file_path  => { strict_type => 1, store => \$map_file_path },
        keep_file_path => { strict_type => 1, store => \$keep_file_path },
        read_freqfile_path =>
          { strict_type => 1, store => \$read_freqfile_path },
        const_fid       => { strict_type => 1, store => \$const_fid },
        sex_check_min_f => { strict_type => 1, store => \$sex_check_min_f },
        vcf_half_call   => {
            allow       => [ undef, qw{ error haploid missing reference } ],
            strict_type => 1,
            store => \$vcf_half_call
        },
        no_fail => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$no_fail
        },
        check_sex => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$check_sex
        },
        y_only => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$y_only
        },
        vcf_require_gt => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$vcf_require_gt
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
        recode => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$recode
        },
        make_just_fam => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$make_just_fam
        },
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
        freqx => {
            default     => 0,
            allow       => [ 0, 1 ],
            strict_type => 1,
            store       => \$freqx
        },
        allow_no_sex => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$allow_no_sex
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = q{plink2};

    ## Options
    if ($binary_fileset_prefix) {

        push @commands, q{--bfile} . $SPACE . $binary_fileset_prefix;
    }

    if ($split_x) {

        push @commands, q{--split-x} . $SPACE . $split_x;

        ## Option modifier
        if ($no_fail) {

            push @commands, q{no-fail};
        }
    }

    if ($make_bed) {

        push @commands, q{--make-bed};
    }

    if ($make_just_fam) {

        push @commands, q{--make-just-fam};
    }

    if ($fam_file_path) {

        push @commands, q{--fam} . $SPACE . $fam_file_path;
    }

    if ($ped_file_path) {

        push @commands, q{--ped} . $SPACE . $ped_file_path;
    }

    if ($map_file_path) {

        push @commands, q{--map} . $SPACE . $map_file_path;
    }

    if ($keep_file_path) {

        push @commands, q{--keep} . $SPACE . $keep_file_path;
    }

    if ($recode) {

        push @commands, q{--recode};
    }

    if ($vcf_require_gt) {

        push @commands, q{--vcf-require-gt};
    }

    if ($vcf_half_call) {

        push @commands, q{--vcf-half-call} . $SPACE . $vcf_half_call;
    }

    if ($set_missing_var_ids) {

        push @commands,
          q{--set-missing-var-ids} . $SPACE . $set_missing_var_ids;
    }

    if ($indep) {

        push @commands, q{--indep};

        ## Option modifier
        if ($indep_window_size) {

            push @commands, $indep_window_size;
        }
        if ($indep_step_size) {

            push @commands, $indep_step_size;
        }
        if ($indep_vif_threshold) {

            push @commands, $indep_vif_threshold;
        }
    }

    if ($check_sex) {

        push @commands, q{--check-sex};

        ## Option modifier
        if ($y_only) {

            push @commands, q{y-only};
        }
        if ($sex_check_min_f) {

            push @commands, $sex_check_min_f;
        }
    }

    if ($het) {

        push @commands, q{--het};

        ## Option modifier
        if ($small_sample) {

            push @commands, q{small-sample};
        }
    }

    if ($freqx) {

        push @commands, q{--freqx};
    }

    if ($allow_no_sex) {

        push @commands, q{--allow-no-sex};
    }

    if ($read_freqfile_path) {

        push @commands, q{--read-freq} . $SPACE . $read_freqfile_path;
    }

    if ($inbreeding_coefficients) {

        push @commands, q{--ibc};
    }

    if ($cluster) {

        push @commands, q{--cluster};
    }

    if ($matrix) {

        push @commands, q{--matrix};
    }

    if ($extract_file) {

        push @commands, q{--extract} . $SPACE . $extract_file;
    }

    #Limit output to regions
    if ( @{$regions_ref} ) {

        push @commands, q{--chr} . $SPACE . join $COMMA, @{$regions_ref};
    }

    if ($vcffile_path) {

        push @commands, q{--vcf} . $SPACE . $vcffile_path;
    }

    if ($const_fid) {

        push @commands, q{--const-fid} . $SPACE . $const_fid;
    }

    ## Outfile
    if ($outfile_prefix) {

        push @commands, q{--out} . $SPACE . $outfile_prefix;
    }

    push @commands,
      unix_standard_streams(
        {
            stdoutfile_path => $stdoutfile_path,
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
