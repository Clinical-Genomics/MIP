package MIP::Recipes::Install::Gatk;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings qw{ FATAL utf8 };
use warnings;

## CPAN
use autodie qw{ :all };

## MIPs lib/
use MIP::Constants qw{ $DOT $LOG_NAME $NEWLINE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.01;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ gatk_download };
}

sub gatk_download {

## Function : Perl wrapper for gatk.tar.bz2 download. Prerequisite for gatk_register and installation into conda environment.
## Returns  : $gatk_tar_path
## Arguments: $gatk_version => Gatk version
##          : $FILEHANDLE   => Filehandle to write to
##          : $quiet        => Be quiet
##          : $verbose      => Set verbosity

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $gatk_version;
    my $FILEHANDLE;
    my $quiet;
    my $verbose;

    ## Default(s)

    my $tmpl = {
        gatk_version => {
            allow       => qr/ \A \d+ \z | \A \d+.\d+ \z /xsm,
            defined     => 1,
            required    => 1,
            store       => \$gatk_version,
            strict_type => 1,
        },
        FILEHANDLE => {
            store => \$FILEHANDLE,
        },
        quiet => {
            allow       => [ undef, 0, 1 ],
            store       => \$quiet,
            strict_type => 1,
        },
        verbose => {
            allow       => [ undef, 0, 1 ],
            store       => \$verbose,
            strict_type => 1,
        },
    };
    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Program::Download::Wget qw{ wget };
    use MIP::Script::Utils qw{ create_temp_dir };

    ## Creating temporary install directory
    say {$FILEHANDLE} q{## Create temporary GATK install directory};
    my $temp_dir_path = create_temp_dir( { FILEHANDLE => $FILEHANDLE } );
    say {$FILEHANDLE} $NEWLINE;

    ## Download
    say {$FILEHANDLE} q{## Download GATK};
    ## NOTE: Since there is no version in the url, we can potentially download
    ## incorrect version of gatk. Gatk-register does not work with the git repo
    ## releases (.zip or .tar.gz)
    my $gatk_url = q{https://software.broadinstitute.org/gatk/download/auth?package=GATK};
    my $gatk_tar_path =
      catfile( $temp_dir_path, q{GenomeAnalysisTK-} . $gatk_version . $DOT . q{tar.bz2} );
    wget(
        {
            FILEHANDLE   => $FILEHANDLE,
            outfile_path => $gatk_tar_path,
            quiet        => $quiet,
            url          => $gatk_url,
            verbose      => $verbose,
        }
    );
    say {$FILEHANDLE} $NEWLINE;
    return $gatk_tar_path;
}

1;
