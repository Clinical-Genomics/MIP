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
use MIP::Constants qw{ $EMPTY_STR $LOG_NAME $PIPE $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ get_binary_version get_executable };
}

sub get_binary_version {

## Function : Get executable/binary versions
## Returns  : %binary_version
## Arguments: $binary_info_href => Binary_Info_Href object

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $binary_info_href;

    my $tmpl = {
        binary_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$binary_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::Executable qw{ get_executable };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my %binary_version;
    my %executable = get_executable( {} );

  BINARY:
    while ( my ( $binary, $binary_path ) = each %{$binary_info_href} ) {

        ## No information on how to get version for this binary - skip
        next BINARY if ( not exists $executable{$binary} );

        ## Unpack
        my $version_cmd    = $executable{$binary}{version_cmd};
        my $version_regexp = $executable{$binary}{version_regexp};

        ## Get
        my @version_cmds = _build_version_cmd(
            {
                binary_path    => $binary_path,
                version_cmd    => $version_cmd,
                version_regexp => $version_regexp,
            }
        );

        ## Call binary and parse output to generate version
        my %cmd_output = system_cmd_call(
            {
                command_string => join $SPACE,
                @version_cmds,
            }
        );
        if ( not $cmd_output{output}[0] ) {

            $log->warn(qq{Could not find version for binary: $binary});
        }
        ## Set binary version
        $binary_version{$binary} = $cmd_output{output}[0];
    }
    return %binary_version;
}

sub get_executable {

## Function : Define the executable features and return them
## Returns  : %{ $executable{$executable_name} } or %executable
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
        mip => {
            version_cmd => q{version},
            version_regexp =>
q?'my ($version) = /mip\s+version\s+(\S+)/xms; if($version) {print $version;last;}'?,
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
        trim_galore => {
            version_cmd => q{--version},
            version_regexp =>
q?'my ($version) = /version\s(\S+)/xms; if($version) {print $version;last;}'?,
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
        just_to_enable_testing => {
            version_cmd    => q{bwa 2>&1 >/dev/null},
            version_regexp => q?'print $version;last;'?,
        },
    );

    if ( defined $executable_name and exists $executable{$executable_name} ) {

        return %{ $executable{$executable_name} };
    }
    return %executable;
}

sub _build_version_cmd {

## Function : Execute mip vercollect to get executable versions
## Returns  : @version_cmds
## Arguments: $binary_path    => Executables (binary) file path
##          : $version_cmd    => Version command line option
##          : $version_regexp => Version reg exp to get version from system call

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $binary_path;
    my $version_cmd;
    my $version_regexp;

    my $tmpl = {
        binary_path => {
            defined     => 1,
            required    => 1,
            store       => \$binary_path,
            strict_type => 1,
        },
        version_cmd => {
            store       => \$version_cmd,
            strict_type => 1,
        },
        version_regexp => {
            defined     => 1,
            required    => 1,
            store       => \$version_regexp,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Language::Perl qw{ perl_nae_oneliners };
    use MIP::Unix::System qw{ system_cmd_call };

    ## Get perl wrapper around regexp
    my @perl_commands = perl_nae_oneliners(
        {
            oneliner_cmd => $version_regexp,
        }
    );

    my @version_cmds = ($binary_path);

    ## Allow for binary without any option to get version
    if ($version_cmd) {

        push @version_cmds, $version_cmd;
    }
    push @version_cmds, ( $PIPE, @perl_commands );
    return @version_cmds;
}

1;
