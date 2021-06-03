package MIP::Environment::Executable;

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
use MIP::Constants qw{ $EMPTY_STR $LOG_NAME $NEWLINE $PIPE $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      build_binary_version_cmd
      get_binary_version_cmd
      get_executable
      get_executable_base_command
      get_binary_version
      write_binaries_versions
    };
}

sub build_binary_version_cmd {

## Function : Build binary version commands
## Returns  : @version_cmds
## Arguments: $binary_path            => Executables (binary) file path
##          : $stdoutfile_path_append => Append stdout info to file path
##          : $use_container          => Use container perl
##          : $version_cmd            => Version command line option
##          : $version_regexp         => Version reg exp to get version from system call

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $binary_path;
    my $stdoutfile_path_append;
    my $version_cmd;
    my $version_regexp;

    ## Default(s)
    my $use_container;

    my $tmpl = {
        binary_path => {
            defined     => 1,
            required    => 1,
            store       => \$binary_path,
            strict_type => 1,
        },
        stdoutfile_path_append => {
            store       => \$stdoutfile_path_append,
            strict_type => 1,
        },
        use_container => => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$use_container,
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

    ## Get perl wrapper around regexp
    my @perl_commands = perl_nae_oneliners(
        {
            oneliner_cmd           => $version_regexp,
            stdoutfile_path_append => $stdoutfile_path_append,
            use_container          => $use_container,
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

sub get_binary_version_cmd {

## Function : Get binary version cmd
## Returns  : @version_cmds
## Arguments: $binary                 => Binary to get version of
##          : $binary_cmd             => Path to binary
##          : $stdoutfile_path_append => Append stdout info to file path
##          : $use_container          => Use container perl

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $binary;
    my $binary_cmd;
    my $stdoutfile_path_append;

    ## Default(s)
    my $use_container;

    my $tmpl = {
        binary => {
            defined     => 1,
            required    => 1,
            store       => \$binary,
            strict_type => 1,
        },
        binary_cmd => {
            defined     => 1,
            required    => 1,
            store       => \$binary_cmd,
            strict_type => 1,
        },
        stdoutfile_path_append => { store => \$stdoutfile_path_append, strict_type => 1, },
        use_container          => => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$use_container,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my %executable = get_executable( { executable_name => $binary, } );

    return if ( not %executable );

    ## Get version command
    my @version_cmds = build_binary_version_cmd(
        {
            binary_path            => $binary_cmd,
            stdoutfile_path_append => $stdoutfile_path_append,
            use_container          => $use_container,
            version_cmd            => $executable{version_cmd},
            version_regexp         => $executable{version_regexp},
        }
    );

    return @version_cmds;
}

sub get_binary_version {

## Function : Get version for executable/binary
## Returns  : $binary_version or %process_return
## Arguments: $binary                   => Binary to get version of
##          : $binary_path              => Path to binary
##          : $capture_version_cmd_href => Capture version cmd hash {REF}
##          : $return_process_hash      => Return full output hash

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $binary;
    my $binary_path;
    my $capture_version_cmd_href;

    ## Default(s)
    my $return_process_hash;

    my $tmpl = {
        binary => {
            defined     => 1,
            required    => 1,
            store       => \$binary,
            strict_type => 1,
        },
        binary_path => {
            defined     => 1,
            required    => 1,
            store       => \$binary_path,
            strict_type => 1,
        },
        capture_version_cmd_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$capture_version_cmd_href,
            strict_type => 1,
        },
        return_process_hash => {
            allow       => [ 0, 1 ],
            default     => 0,
            store       => \$return_process_hash,
            strict_type => => 1,
        }
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Environment::Child_process qw{ child_process };

    ## Get version command
    my @version_cmds = build_binary_version_cmd(
        {
            binary_path    => $binary_path,
            version_cmd    => $capture_version_cmd_href->{version_cmd},
            version_regexp => $capture_version_cmd_href->{version_regexp},
        }
    );

    ## Call binary and parse output to generate version
    my %process_return = child_process(
        {
            commands_ref => [ join $SPACE, @version_cmds ],
            process_type => q{ipc_cmd_run},
        }
    );

    return %process_return if ($return_process_hash);

    return $process_return{stdouts_ref}[0];
}

sub get_executable_base_command {

## Function : Get executable base command with or without container manager commands
## Returns  : $base_command
## Arguments: $base_command => Executable standard base command

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $base_command;

    my $tmpl = {
        base_command => {
            defined     => 1,
            required    => 1,
            store       => \$base_command,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Constants qw{ %CONTAINER_CMD };

    return exists $CONTAINER_CMD{$base_command}
      ? $CONTAINER_CMD{$base_command}
      : $base_command;
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
        arriba => {
            version_cmd    => q{-h},
            version_regexp =>
              q?'my ($version) = /\AVersion:\s(\S+)/xms; if($version) {print $version;last;}'?,
        },
        q{bam2wig.py} => {
            version_cmd    => q{--version},
            version_regexp =>
              q?'my ($version) = /bam2wig.py\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        q{bam_stat.py} => {
            version_cmd    => q{--version},
            version_regexp =>
              q?'my ($version) = /bam_stat.py\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        bcftools => {
            version_cmd    => q{2>&1 >/dev/null},
            version_regexp =>
q?'my ($version) = /Version:\s+(.*)/xms; if($version) {chomp $version;print $version;last;}'?,
        },
        bedtools => {
            version_cmd    => q{--version},
            version_regexp =>
              q?'my ($version) = /bedtools\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        bgzip => {
            version_cmd    => q{-h 2>&1 },
            version_regexp =>
              q?'my ($version) = /Version:\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        bwa => {
            version_cmd    => q{2>&1 >/dev/null},
            version_regexp =>
              q?'my ($version) = /Version:\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        q{bwa-mem2} => {
            version_cmd    => q{ version 2>&1 >/dev/null},
            version_regexp =>
              q?'my ($version) = /Version:\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        chanjo => {
            version_cmd    => q{--version},
            version_regexp =>
              q?'my ($version) = /version\s(\S+)/xms; if($version) {print $version;last;}'?,
        },
        q{configManta.py} => {
            version_cmd    => q{--version},
            version_regexp => q?'chomp;print $_;last;'?,
        },
        fastqc => {
            version_cmd    => q{--version},
            version_regexp =>
              q?'my ($version) = /FastQC\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        ExpansionHunter => {
            version_cmd    => q{--version 2>&1 >/dev/null},
            version_regexp =>
q?'my ($version) = /Hunter\s+(v\d+.\d+.\d+)/xms; if($version) {print $version;last;}'?,
        },
        gatk => {
            version_cmd    => q{--java-options "-Xmx1G" --version},
            version_regexp =>
              q?'my ($version) = /\(GATK\)\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        q{geneBody_coverage2.py} => {
            version_cmd    => q{--version},
            version_regexp =>
q?'my ($version) = /geneBody_coverage2.py\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        genmod => {
            version_cmd    => q{--version},
            version_regexp =>
q?'my ($version) = /genmod\s+version:\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        gffcompare => {
            version_cmd    => q{--version 2>&1 >/dev/null},
            version_regexp =>
              q?'my ($version) = /gffcompare\sv(\S+)/xms; if($version) {print $version;last;}'?,
        },
        q{grep} => {
            version_cmd    => q{--version},
            version_regexp =>
              q?'my ($version) = /\(GNU\s+grep\)\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        gzip => {
            version_cmd    => q{--version},
            version_regexp =>
              q?'my ($version) = /gzip\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        q{infer_exeperiment.py} => {
            version_cmd    => q{--version},
            version_regexp =>
q?'my ($version) = /infer_exeperiment.py\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        q{inner_distance.py} => {
            version_cmd    => q{--version},
            version_regexp =>
q?'my ($version) = /inner_distance.py\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        q{junction_annotation.py} => {
            version_cmd    => q{--version},
            version_regexp =>
q?'my ($version) = /junction_annotation.py\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        mip => {
            version_cmd    => q{version},
            version_regexp =>
              q?'my ($version) = /mip\s+version\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        multiqc => {
            version_cmd    => q{--version},
            version_regexp =>
              q?'my ($version) = /version\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        peddy => {
            version_cmd    => q{--version},
            version_regexp =>
              q?'my ($version) = /version\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        preseq => {
            version_cmd    => q{2>&1 >/dev/null},
            version_regexp =>
              q?'my ($version) = /Version:\s+(\S+)/xms; if ($version) {print $version; last;}'?,
        },
        picard => {
            version_cmd    => q{java -jar /usr/picard/picard.jar BamIndexStats 2>&1 >/dev/null},
            version_regexp =>
              q?'my ($version) = /Version:(\S+)/xms; if($version) {print $version;last;}'?,
        },
        pigz => {
            version_cmd    => q{--version 2>&1 >/dev/null},
            version_regexp =>
              q?'my ($version) = /pigz\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        plink2 => {
            version_cmd    => q{--version},
            version_regexp =>
q?'my ($version) = /PLINK\s+(.*)/xms; if($version) {chomp $version;print $version;last;}'?,
        },
        q{read_distribution.py} => {
            version_cmd    => q{--version},
            version_regexp =>
q?'my ($version) = /read_distribution.py\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        q{read_duplication.py} => {
            version_cmd    => q{--version},
            version_regexp =>
q?'my ($version) = /read_duplication.py\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        rhocall => {
            version_cmd    => q{--version},
            version_regexp =>
              q?'my ($version) = /version\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        salmon => {
            version_cmd    => q{--version},
            version_regexp =>
              q?'my ($version) = /salmon\s+(\S+)/xms; if ($version) {print $version; last;}'?,
        },
        sambamba => {
            version_cmd    => q{--version 2>&1 >/dev/null},
            version_regexp =>
              q?'my ($version) = /sambamba\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        samtools => {
            version_cmd    => q{2>&1 >/dev/null},
            version_regexp =>
q?'my ($version) = /Version:\s+(.*)/xms; if($version) {chomp $version;print $version;last;}'?,
        },
        sed => {
            version_cmd    => q{--version},
            version_regexp =>
              q?'my ($version) = /\(GNU\s+sed\)\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        STAR => {
            version_cmd    => q{--version},
            version_regexp => q?'my ($version) = /(\S+)/xms; if($version) {print $version; last;}'?,
        },
        q{STAR-Fusion} => {
            version_cmd    => q{--version},
            version_regexp =>
              q?'my ($version) = /version:\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        stranger => {
            version_cmd    => q{--version},
            version_regexp => q?'my ($version) = /(\S+)/xms; if($version) {print $version;last;}'?,
        },
        stringtie => {
            version_cmd    => q{--version},
            version_regexp => q?'my ($version) = /(\S+)/xms; if($version) {print $version;last;}'?
        },
        svdb => {
            version_cmd    => q{--version},
            version_regexp =>
              q?'my ($version) = /SVDB-(\S+)/xms; if($version) {print $version;last;}'?,
        },
        tabix => {
            version_cmd    => q{--version},
            version_regexp =>
              q?'my ($version) = /\(htslib\)\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        telomerecat => {
            version_regexp =>
q?'BEGIN {my $match = 0}; if ($match) {my ($version) = /(\d+.\d+.\d+)/; print $version; last;} $match=1 if (/\A Version: /xms);'?,
        },
        q{TIDDIT.py} => {
            version_cmd    => q{--version},
            version_regexp =>
              q?'my ($version) = /TIDDIT-(\S+)/xms; if($version) {print $version;last;}'?,
        },
        trim_galore => {
            version_cmd    => q{--version},
            version_regexp =>
              q?'my ($version) = /version\s(\S+)/xms; if($version) {print $version;last;}'?,
        },
        upd => {
            version_cmd    => q{--version},
            version_regexp => q?'my ($version) = /(\S+)/xms; if($version) {print $version;last;}'?,
        },
        vcfanno => {
            version_cmd    => q{2>&1 >/dev/null},
            version_regexp =>
              q?'my ($version) = /version\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        vcf2cytosure => {
            version_cmd    => q{-V 2>&1},
            version_regexp =>
              q?'my ($version) = /to\s+cytosure\s+(\S+)/xms; if($version) {print $version;last;}'?,
        },
        varg => {
            version_cmd    => q{--version},
            version_regexp =>
              q?'my ($version) = /version\s(\S+)/xms; if($version) {print $version;last;}'?
        },
        vep => {
            version_regexp =>
              q?'my ($version) = /ensembl-vep\s+:\s(\d+)/xms; if($version) {print $version;last;}'?,
        },
        q{wigToBigWig} => {
            version_cmd    => q{2>&1 >/dev/null},
            version_regexp =>
              q?'my ($version) = /wigToBigWig\sv\s(\S+)/xms; if($version) {print $version;last;}'?,
        },
        just_to_enable_testing => {
            version_cmd    => q{bwa 2>&1 >/dev/null},
            version_regexp => q?'print $version;last;'?,
        },
    );

    if ( defined $executable_name and exists $executable{$executable_name} ) {

        return %{ $executable{$executable_name} };
    }

    ## Missing information on executable
    return if ( defined $executable_name );

    return %executable;
}

sub write_binaries_versions {

## Function : Write executables/binaries versions
## Returns  :
## Arguments: $binary_info_href => Binary info href object
##          : $filehandle       => Filehandle to write to
##          : $outfile_path     => Outfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $binary_info_href;
    my $filehandle;
    my $outfile_path;

    my $tmpl = {
        binary_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$binary_info_href,
            strict_type => 1,
        },
        filehandle => {
            required => 1,
            store    => \$filehandle,
        },
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Language::Perl qw{ perl_nae_oneliners };

  BINARY:
    while ( my ( $binary, $binary_cmd ) = each %{$binary_info_href} ) {

        my @version_cmds = get_binary_version_cmd(
            {
                binary        => $binary,
                binary_cmd    => $binary_cmd,
                use_container => 1,
            }
        );

        ## No information on how to get version for this binary - skip
        next BINARY if ( not @version_cmds );

        my @add_binary_cmds = perl_nae_oneliners(
            {
                oneliner_name          => q{add_binary},
                oneliner_parameter     => $binary,
                stdoutfile_path_append => $outfile_path,
                use_container          => 1,
            }
        );
        push @version_cmds, ( $PIPE, @add_binary_cmds );
        say {$filehandle} join( $SPACE, @version_cmds ), $NEWLINE;
    }
    return;
}

1;
