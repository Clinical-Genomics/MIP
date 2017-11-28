package MIP::Recipes::Install::Gatk;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use File::Spec::Functions qw{ catfile };

## CPANM
use Readonly;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ gatk_download };
}

## Constants
Readonly my $DOT     => q{.};
Readonly my $NEWLINE => qq{\n};

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
    my $quiet;
    my $verbose;
    my $FILEHANDLE;

    ## Default(s)

    my $tmpl = {
        gatk_version => {
            required    => 1,
            defined     => 1,
            allow       => qr/ ^\d+$ | ^\d+.\d+$ /xsm,
            strict_type => 1,
            store       => \$gatk_version
        },
        FILEHANDLE => {
            store => \$FILEHANDLE
        },
        quiet => {
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$quiet
        },
        verbose => {
            allow       => [ undef, 0, 1 ],
            strict_type => 1,
            store       => \$verbose
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
    my $gatk_url =
      q{https://software.broadinstitute.org/gatk/download/auth?package=GATK};
    my $gatk_tar_path =
      catfile( $temp_dir_path,
        q{GenomeAnalysisTK-} . $gatk_version . $DOT . q{tar.bz2} );
    wget(
        {
            url          => $gatk_url,
            FILEHANDLE   => $FILEHANDLE,
            quiet        => $quiet,
            verbose      => $verbose,
            outfile_path => $gatk_tar_path,
        }
    );
    say {$FILEHANDLE} $NEWLINE;
    return $gatk_tar_path;
}

1;
