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
    our $VERSION = 1.00;

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
            allow       => qr{ \A\d+\z }sxm,
            default     => 0,
            store       => \$autosplit,
            strict_type => 1,
        },
        command_line => {
            allow       => qr{ \A\d+\z }sxm,
            default     => 0,
            store       => \$command_line,
            strict_type => 1,
        },
        n => {
            allow       => qr{ \A\d+\z }sxm,
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
## Arguments: $autosplit                 => Turns on autosplit mode when used with a -n or -p
##          : $command_line              => Enter one line of program
##          : $FILEHANDLE                => Filehandle to write to
##          : $n                         => Iterate over filename arguments
##          : $stderrfile_path           => Stderrfile path
##          : $stderrfile_path_append    => Append stderr info to file path
##          : $stdinfile_path            => Stdinfile path
##          : $stdoutfile_path           => Stdoutfile path
##          : $synonyms_grch37_to_grch38 => Predefined oneliner
##          : $synonyms_grch38_to_grch37 => Predefined oneliner

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $FILEHANDLE;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdinfile_path;
    my $stdoutfile_path;

    ## Default(s)
    my $autosplit;
    my $command_line;
    my $n;
    my $synonyms_grch37_to_grch38;
    my $synonyms_grch38_to_grch37;

    my $tmpl = {
        autosplit => {
            allow       => qr{ \A\d+\z }sxm,
            default     => 1,
            store       => \$autosplit,
            strict_type => 1,
        },
        command_line => {
            allow       => qr{ \A\d+\z }sxm,
            default     => 1,
            store       => \$command_line,
            strict_type => 1,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        n => {
            allow       => qr{ \A\d+\z }sxm,
            default     => 1,
            store       => \$n,
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
        synonyms_grch37_to_grch38 => {
            allow       => qr{ \A\d+\z }sxm,
            default     => 0,
            store       => \$synonyms_grch37_to_grch38,
            strict_type => 1,
        },
        synonyms_grch38_to_grch37 => {
            allow       => qr{ \A\d+\z }sxm,
            default     => 0,
            store       => \$synonyms_grch38_to_grch37,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Stores commands depending on input parameters
    my @commands = perl_base(
        {
            autosplit    => $autosplit,
            command_line => $command_line,
            n            => $n,
        }
    );

    ## \' is for xargs
    if ($synonyms_grch37_to_grch38) {

        push @commands,
          q?\'if($_=~s/^M/chrMT/g) {} elsif ($_=~s/^(.+)/chr$1/g) {} print $_\'?;
    }

    ## \' is for xargs
    if ($synonyms_grch38_to_grch37) {

        push @commands,
          q?\'if($_=~s/^chrMT/M/g) {} elsif ($_=~s/^chr(.+)/$1/g) {} print $_\'?;
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
            FILEHANDLE   => $FILEHANDLE,
            separator    => $SPACE,

        }
    );
    return @commands;
}

1;
