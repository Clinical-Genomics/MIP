package MIP::Recipes::Install::Star_fusion;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catfile devnull };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings qw{ FATAL utf8 };
use warnings;

## CPAN
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $AT $DOLLAR_SIGN $DOUBLE_QUOTE $NEWLINE $SINGLE_QUOTE $SPACE };
use MIP::Gnu::Coreutils qw{ gnu_chmod gnu_echo };
use MIP::Language::Shell qw{ build_shebang };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.04;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ setup_star_fusion };
}

sub setup_star_fusion {

## Function : Put executables from STAR-Fusion and its utils in condas path
## Returns  :
## Arguments: $conda_prefix_path    => Conda prefix path
##          : $FILEHANDLE           => Filehandle to write to
##          : $star_fusion_dir_path => Path to STAR-Fusion directory

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $conda_prefix_path;
    my $FILEHANDLE;
    my $star_fusion_dir_path;

    my $tmpl = {
        conda_prefix_path => {
            defined     => 1,
            required    => 1,
            store       => \$conda_prefix_path,
            strict_type => 1,
        },
        FILEHANDLE => {
            defined  => 1,
            required => 1,
            store    => \$FILEHANDLE,
        },
        star_fusion_dir_path => {
            required    => 1,
            store       => \$star_fusion_dir_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @shebang = build_shebang(
        {
            bash_bin_path => catfile( dirname( dirname( devnull() ) ), qw{ bin bash } ),
        }
    );

    my @util_executables = qw{
      gtf_file_to_feature_seqs.pl
      make_super_locus.pl
      remove_long_intron_readthru_transcripts.pl
      restrict_genome_to_chr_entries.pl
    };

    my %executable = map {
        $_ => catfile( $star_fusion_dir_path, qw {ctat-genome-lib-builder util }, $_ )
    } @util_executables;

    ## The addition of prep_genome_lib.pl to bin via conda seems to be broken in STAR-Fusion version 1.7.0-1.
    $executable{q{prep_genome_lib.pl}} =
      catfile( $star_fusion_dir_path, qw{ ctat-genome-lib-builder prep_genome_lib.pl } );

    say {$FILEHANDLE} q{## Creating proxy executables};
  MOCK_EXECUTABLE:
    foreach my $proxy_executable ( keys %executable ) {

        my $file_content =
            $executable{$proxy_executable}
          . $SPACE
          . $DOUBLE_QUOTE
          . $DOLLAR_SIGN
          . $AT
          . $DOUBLE_QUOTE;
        my $outfile_path = catfile( $conda_prefix_path, q{bin}, $proxy_executable );
        gnu_echo(
            {
                FILEHANDLE     => $FILEHANDLE,
                outfile_path   => $outfile_path,
                string_wrapper => $SINGLE_QUOTE,
                strings_ref    => \@shebang,
            }
        );
        print {$FILEHANDLE} $NEWLINE;
        gnu_echo(
            {
                FILEHANDLE             => $FILEHANDLE,
                no_trailing_newline    => 1,
                stdoutfile_path_append => $outfile_path,
                string_wrapper         => $SINGLE_QUOTE,
                strings_ref            => [$file_content],
            }
        );
        print {$FILEHANDLE} $NEWLINE;
        gnu_chmod(
            {
                file_path  => $outfile_path,
                FILEHANDLE => $FILEHANDLE,
                permission => q{a+x},
            }
        );
        print {$FILEHANDLE} $NEWLINE;
    }
    return;
}

1;

