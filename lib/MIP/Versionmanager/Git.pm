package MIP::Versionmanager::Git;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use FindBin qw{ $Bin };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ git_clone git_checkout};
}

## Constants
Readonly my $SPACE => q{ };

sub git_clone {

## Function : Perl wrapper for git clone. Based on version 2.13.3.
## Returns  : @commands
## Arguments: $url                    => Url to use for clone
##          : $outdir_path            => Path to output directory, must be empty
##          : $stdoutfile_path        => Stdoutfile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $filehandle             => Filehandle to write to
##          : $verbose                => Verbose output
##          : $quiet                  => Less chatty output

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $url;
    my $outdir_path;
    my $stdoutfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $filehandle;

    ## Default(s)
    my $verbose;
    my $quiet;

    my $tmpl = {
        url => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$url
        },
        outdir_path => {
            strict_type => 1,
            store       => \$outdir_path,
        },
        stdoutfile_path => {
            strict_type => 1,
            store       => \$stdoutfile_path
        },
        stderrfile_path => {
            strict_type => 1,
            store       => \$stderrfile_path
        },
        stderrfile_path_append => {
            strict_type => 1,
            store       => \$stderrfile_path_append
        },
        filehandle => {
            store => \$filehandle
        },
        verbose => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$verbose,
        },
        quiet => {
            default     => 0,
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$quiet,
        },

    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ git clone };

    if ($verbose) {
        push @commands, q{--verbose};
    }

    if ($quiet) {
        push @commands, q{--quiet};
    }

    push @commands, $url;

    if ($outdir_path) {
        push @commands, $outdir_path;
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
            filehandle   => $filehandle,
        }
    );

    return @commands;
}

sub git_checkout {

## Function : Perl wrapper for git checkout. Based on version 2.13.3.
## Returns  : @commands
## Arguments: $branch                 => Branch to checkout
##          : $stdoutfile_path        => Stdoutfile path
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $filehandle             => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $branch;
    my $stdoutfile_path;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $filehandle;

    ## Default(s)

    my $tmpl = {
        branch =>
          { required => 1, defined => 1, strict_type => 1, store => \$branch },
        stdoutfile_path => { strict_type => 1, store => \$stdoutfile_path },
        stderrfile_path => { strict_type => 1, store => \$stderrfile_path },
        stderrfile_path_append =>
          { strict_type => 1, store => \$stderrfile_path_append },
        filehandle => { store => \$filehandle },

    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Stores commands depending on input parameters
    my @commands = qw{ git checkout };

    push @commands, $branch;

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
            filehandle   => $filehandle,
        }
    );

    return @commands;
}

1;
