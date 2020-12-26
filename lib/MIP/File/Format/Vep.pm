package MIP::File::Format::Vep;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Path qw{ make_path };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $LOG_NAME $NEWLINE %PRIMARY_CONTIG $SPACE $TAB };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ create_vep_synonyms_file };
}

sub create_vep_synonyms_file {

## Function : Create the synonyms file for VEP option '--synonyms'
## Returns  : undef or $outfile_path
## Arguments: $outfile_path => Outfile path to write to
##          : $version      => Human genome version {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $outfile_path;
    my $version;

    ## Default(s)

    my $tmpl = {
        outfile_path => {
            defined     => 1,
            required    => 1,
            store       => \$outfile_path,
            strict_type => 1,
        },
        version => {
            defined     => 1,
            required    => 1,
            store       => \$version,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Contigs qw{ get_contig_set };

    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my %contig_synonyms_map = get_contig_set(
        {
            contig_set => q{synonyms_map},
            version    => $version,
        }
    );

    return if ( not %contig_synonyms_map );

    my @contigs = get_contig_set(
        {
            contig_set => q{contigs},
            version    => $version,
        }
    );

    $log->info( q{Creating VEP synonyms file: } . $outfile_path, $NEWLINE );

    ## Create dir if it does not exists
    make_path( dirname($outfile_path) );

    open my $filehandle, q{>}, $outfile_path
      or $log->logdie(qq{Cannot open $outfile_path: $ERRNO });

  CONTIG:
    foreach my $primary_contig (@contigs) {

        my $synonymous_contig = $contig_synonyms_map{$primary_contig};

        ## No defined synonym
        next CONTIG if ( not $synonymous_contig );

        say {$filehandle} $primary_contig . $TAB . $synonymous_contig;
    }
    $log->info( q{Wrote: } . $outfile_path, $NEWLINE );
    close $filehandle;
    return $outfile_path;
}

1;
