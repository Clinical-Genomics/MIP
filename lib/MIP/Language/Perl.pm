package MIP::Language::Perl;

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
use MIP::Constants qw{ $DASH $SPACE };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.03;

    our @EXPORT_OK = qw{ perl_base perl_nae_oneliners };
}

sub perl_base {

## Function : Perl base and switches
## Returns  :
## Arguments: $autosplit    => Turns on autosplit mode when used with a -n or -p
##          : $command_line => Enter one line of program
##          : $n            => Iterate over filename arguments

    my ($arg_href) = @_;

    ## Flatten argument(s)

    ## Default(s)
    my $autosplit;
    my $command_line;
    my $n;

    my $tmpl = {
        autosplit => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$autosplit,
            strict_type => 1,
        },
        command_line => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$command_line,
            strict_type => 1,
        },
        n => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$n,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = qw{ perl };

    if ($n) {

        push @commands, q{-n};
    }
    if ($autosplit) {

        push @commands, q{-a};
    }
    if ($command_line) {

        push @commands, q{-e};
    }

    return @commands;
}

sub perl_nae_oneliners {

## Function : Return predifined one liners
## Returns  : @commands
## Arguments: $autosplit              => Turns on autosplit mode when used with a -n or -p
##          : $command_line           => Enter one line of program
##          : $filehandle             => Filehandle to write to
##          : $n                      => Iterate over filename arguments
##          : $oneliner_cmd           => Command to execute
##          : $oneliner_name          => Perl oneliner name
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdinfile_path         => Stdinfile path
##          : $stdoutfile_path        => Stdoutfile path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $oneliner_cmd;
    my $oneliner_name;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdinfile_path;
    my $stdoutfile_path;

    ## Default(s)
    my $autosplit;
    my $command_line;
    my $n;

    my $tmpl = {
        autosplit => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$autosplit,
            strict_type => 1,
        },
        command_line => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$command_line,
            strict_type => 1,
        },
        filehandle => {
            store => \$filehandle,
        },
        n => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$n,
            strict_type => 1,
        },
        oneliner_cmd => {
            store       => \$oneliner_cmd,
            strict_type => 1,
        },
        oneliner_name => {
            store       => \$oneliner_name,
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
        stdinfile_path  => { store => \$stdinfile_path, strict_type => 1, },
        stdoutfile_path => {
            store       => \$stdoutfile_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Oneliner dispatch table
    my %oneliner = (
        write_contigs_size_file   => \&_write_contigs_size_file,
        get_dict_contigs          => \&_get_dict_contigs,
        get_vep_version           => \&_get_vep_version,
        synonyms_grch37_to_grch38 => \&_synonyms_grch37_to_grch38,
        synonyms_grch38_to_grch37 => \&_synonyms_grch38_to_grch37,
    );

    ## Stores commands depending on input parameters
    my @commands = perl_base(
        {
            autosplit    => $autosplit,
            command_line => $command_line,
            n            => $n,
        }
    );

    if ( defined $oneliner_name and exists $oneliner{$oneliner_name} ) {

        push @commands, $oneliner{$oneliner_name}->();
    }
    elsif ($oneliner_cmd) {

        push @commands, $oneliner_cmd;
    }

    push @commands,
      unix_standard_streams(
        {
            stderrfile_path        => $stderrfile_path,
            stderrfile_path_append => $stderrfile_path_append,
            stdinfile_path         => $stdinfile_path,
            stdoutfile_path        => $stdoutfile_path,
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

sub _get_vep_version {

## Function : Return predifined one liners for getting vep version
## Returns  : $get_vep_version
## Arguments:

    my ($arg_href) = @_;

    my $get_vep_version = q?'if($_=~/ensembl-vep\s+:\s(\d+)/xms) {print $1;}'?;

    return $get_vep_version;
}

sub _get_dict_contigs {

## Function : Return predifined one liners to get contig names from dict file
## Returns  : $get_dict_contigs
## Arguments:

    my ($arg_href) = @_;

    # Find contig line
    my $get_dict_contigs = q?'if($F[0]=~/^\@SQ/) { ?;

    # Collect contig name
    $get_dict_contigs .= q? if($F[1]=~/SN\:(\S+)/) { ?;

    # Alias capture
    $get_dict_contigs .= q?my $contig_name = $1; ?;

    # Write to STDOUT
    $get_dict_contigs .= q?print $contig_name, q{,};} }'?;

    return $get_dict_contigs;
}

sub _synonyms_grch37_to_grch38 {

## Function : Return predifined one liners to modify chr prefix from genome version 37 to 38
## Returns  : $modify_chr_prefix
## Arguments:

    my ($arg_href) = @_;

    ## Add "chr" prefix to chromosome name and rename "M" to "MT"
    my $modify_chr_prefix = q?\'if($_=~s/^M/chrMT/g) {} ?;

    ## Add "chr" prefix to chromosome name
    $modify_chr_prefix .= q?elsif ($_=~s/^(.+)/chr$1/g) {} ?;

## Print line
    $modify_chr_prefix .= q?print $_\'?;

    return $modify_chr_prefix;
}

sub _synonyms_grch38_to_grch37 {

## Function : Return predifined one liners to modify chr prefix from genome version 37 to 38
## Returns  : $modify_chr_prefix
## Arguments:

    my ($arg_href) = @_;

## Remove "chr" prefix from chromosome name and rename "MT" to "M"
    my $modify_chr_prefix = q?\'if($_=~s/^chrMT/M/g) {} ?;

## Remove "chr" prefix from chromosome name
    $modify_chr_prefix .= q?elsif ($_=~s/^chr(.+)/$1/g) {} ?;

## Print line
    $modify_chr_prefix .= q?print $_\'?;

    return $modify_chr_prefix;
}

sub _write_contigs_size_file {

## Function : Return predifined one liners to write contig names and length from fai file
## Returns  : $write_contigs_size
## Arguments:

    my ($arg_href) = @_;

    ## Contig name ($F[0]), contig length ($F[1])
    my $write_contigs_size = q?'say STDOUT $F[0] . "\t" . $F[1] '?;

    return $write_contigs_size;
}

1;
