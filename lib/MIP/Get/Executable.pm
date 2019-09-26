package MIP::Get::Executable;

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
use MIP::Constants qw{ $EMPTY_STR $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ get_executable };
}

sub get_executable {

## Function :
## Returns  :
## Arguments: $executable_name => Executable name

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $executable_name;

    my $tmpl = {
        executable_name => {
            store       => \$executable_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %executable = (
        bcftools => {
            version_cmd => q{2>&1 >/dev/null},
            version_regexp =>
q?'my ($version) = /Version:\s+(.*)/xms; if($version) {chomp $version;print $version;last;}'?,
        },
        bedtools => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /bedtools\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        bgzip => {
            version_cmd => q{-h 2>&1 >/dev/null},
            version_regexp =>
q?'my ($version) = /Version:\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        bwa => {
            version_cmd => q{2>&1 >/dev/null},
            version_regexp =>
q?'my ($version) = /Version:\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        chanjo => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /version\s(\S+)/xms; if($version) {print $version;last;}'?,
        },
        q{configManta.py} => {
            version_cmd    => q{--version},
            version_regexp => q?'chomp;print $_;last;'?,
        },
        fastqc => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /FastQC\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        ExpansionHunter => {
            version_cmd => q{--version 2>&1 >/dev/null},
            version_regexp =>
q?'my ($version) = /Hunter\s+(v\d+.\d+.\d+)/xms; if($version) {print $version;last;}'?,
        },
        gatk => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /\(GATK\)\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        genmod => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /genmod\s+version:\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        q{grep} => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /\(GNU\s+grep\)\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        gzip => {
            version_cmd => q{--version},
            version_regexp =>
              q?'my ($version) = /gzip\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        multiqc => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /version\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        peddy => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /version\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        picard => {
            version_cmd => q{BamIndexStats 2>&1 >/dev/null},
            version_regexp =>
q?'my ($version) = /Version:\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        plink2 => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /PLINK\s+(.*)/xms; if($version) {chomp $version;print $version;last;}'?,
        },
        rhocall => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /version\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        sambamba => {
            version_cmd => q{--version 2>&1 >/dev/null},
            version_regexp =>
q?'my ($version) = /sambamba\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        samtools => {
            version_cmd => q{2>&1 >/dev/null},
            version_regexp =>
q?'my ($version) = /Version:\s+(.*)/xms; if($version) {chomp $version;print $version;last;}'?,
        },
        sed => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /\(GNU\s+sed\)\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        stranger => {
            version_cmd => q{--version},
            version_regexp =>
              q?'my ($version) = /(\S+)/xms; if($version) {print $version;last;}'?,
        },
        svdb => {
            version_cmd => q{--version},
            version_regexp =>
              q?'my ($version) = /SVDB-(\S+)/xms; if($version) {print $version;last;}'?,
        },
        tabix => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /\(htslib\)\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        q{TIDDIT.py} => {
            version_cmd => q{--version},
            version_regexp =>
              q?'my ($version) = /TIDDIT-(\S+)/xms; if($version) {print $version;last;}'?,
        },
        upd => {
            version_cmd => q{--version},
            version_regexp =>
              q?'my ($version) = /(\S+)/xms; if($version) {print $version;last;}'?,
        },
        variant_integrity => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /version\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        vcfanno => {
            version_cmd => q{2>&1 >/dev/null},
            version_regexp =>
q?'my ($version) = /version\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        vcf2cytosure => {
            version_cmd => q{-V 2>&1 >/dev/null},
            version_regexp =>
q?'my ($version) = /to\s+cytosure\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        vep => {
            version_regexp =>
q?'my ($version) = /ensembl-vep\s+:\s(\d+)/xms; if($version) {print $version;last;}'?,
        },
        vt => {
            version_cmd => q{normalize 2>&1 >/dev/null},
            version_regexp =>
q?'my ($version) = /normalize\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
    );

    if ( defined $executable_name and exists $executable{$executable_name} ) {

        return %{ $executable{$executable_name} };
    }
    return %executable;
}

1;
