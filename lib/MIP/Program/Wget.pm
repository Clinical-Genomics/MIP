package MIP::Program::Wget;

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
use MIP::Constants qw{ $EQUALS $SPACE };
use MIP::Unix::Standard_streams qw{ unix_standard_streams };
use MIP::Unix::Write_to_file qw{ unix_write_to_file };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ wget };
}

sub wget {

## Function : Perl wrapper for writing wget recipe to $filehandle or return commands array. Based on GNU Wget 1.12, a non-interactive network retriever.
## Returns  : @commands
## Arguments: $continue               => Resume getting a partially-downloaded file
##          : $filehandle             => Filehandle to write to
##          : $outfile_path           => Outfile path. Write documents to FILE
##          : $quiet                  => Quiet (no output)
##          : $read_timeout           => Set the read timeout to SECS
##          : $retry_connrefused      => Retry even if connection is refused
##          : $stderrfile_path        => Stderrfile path
##          : $stderrfile_path_append => Append stderr info to file path
##          : $stdoutfile_path        => Stdoutfile path
##          : $timeout                => Set all timeout values to SECONDS
##          : $tries                  => Set number of retries to NUMBER (0 unlimits)
##          : $url                    => Url to use for download
##          : $user                   => User name
##          : $verbose                => Verbosity
##          : $wait_retry             => wait 1..SECONDS between retries of a retrieval

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $filehandle;
    my $outfile_path;
    my $read_timeout;
    my $stderrfile_path;
    my $stderrfile_path_append;
    my $stdoutfile_path;
    my $timeout;
    my $tries;
    my $url;
    my $user;
    my $wait_retry;

    ## Default(s)
    my $continue;
    my $quiet;
    my $retry_connrefused;
    my $verbose;

    my $tmpl = {
        continue => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$continue,
            strict_type => 1,
        },
        filehandle   => { store => \$filehandle, },
        outfile_path => { store => \$outfile_path, strict_type => 1, },
        quiet        => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$quiet,
            strict_type => 1,
        },
        read_timeout => {
            allow       => [ undef, qr{ \A\d+\z }sxm ],
            store       => \$read_timeout,
            strict_type => 1,
        },
        retry_connrefused => {
            allow       => [ undef, 0, 1 ],
            default     => 0,
            store       => \$retry_connrefused,
            strict_type => 1,
        },
        stdoutfile_path        => { store => \$stdoutfile_path,        strict_type => 1, },
        stderrfile_path        => { store => \$stderrfile_path,        strict_type => 1, },
        stderrfile_path_append => { store => \$stderrfile_path_append, strict_type => 1, },
        timeout                => {
            allow       => [ undef, qr{ \A\d+\z }sxm ],
            store       => \$timeout,
            strict_type => 1,
        },
        tries => {
            allow       => [ undef, qr{ \A\d+\z }sxm ],
            store       => \$tries,
            strict_type => 1,
        },
        url  => { defined => 1, required => 1, store => \$url, strict_type => 1, },
        user => {
            store       => \$user,
            strict_type => 1,
        },
        verbose => {
            allow       => [ undef, 0, 1 ],
            default     => 1,
            store       => \$verbose,
            strict_type => 1,
        },
        wait_retry => {
            allow       => [ undef, qr{ \A\d+\z }sxm ],
            store       => \$wait_retry,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @commands = qw{ wget };

    if ($quiet) {

        push @commands, q{--quiet};
    }
    if ($verbose) {

        push @commands, q{--verbose};
    }
    else {

        push @commands, q{--no-verbose};
    }

    if ($continue) {

        push @commands, q{--continue};
    }
    if ($retry_connrefused) {

        push @commands, q{--retry-connrefused};
    }
    if ($wait_retry) {

        push @commands, q{--waitretry} . $EQUALS . $wait_retry;
    }

    if ($read_timeout) {

        push @commands, q{--read-timeout} . $EQUALS . $read_timeout;
    }

    if ($timeout) {

        push @commands, q{--timeout} . $EQUALS . $timeout;
    }

    if ($tries) {

        push @commands, q{--tries} . $EQUALS . $tries;
    }
    if ($user) {

        push @commands, q{--user} . $EQUALS . $user;
    }

    ## URL
    push @commands, $url;

    ## Outfile
    if ($outfile_path) {

        push @commands, q{-O} . $SPACE . $outfile_path;
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
            filehandle   => $filehandle,
            separator    => $SPACE,
        }
    );

    return @commands;
}

1;
