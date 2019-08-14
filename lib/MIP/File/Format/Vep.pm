package MIP::File::Format::Vep;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Path qw{ make_path };
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
use MIP::Constants qw{ $NEWLINE %PRIMARY_CONTIG $SPACE $TAB };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ create_vep_synonyms_file };
}

sub create_vep_synonyms_file {

## Function : Create the synonyms file for VEP option '--synonyms'
## Returns  : undef or $outfile_path
## Arguments: $log          => Log object
##          : $outfile_path => Outfile path to write to
##          : $version      => Human genome version {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $log;
    my $outfile_path;
    my $version;

    ## Default(s)

    my $tmpl = {
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
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

    ## Make local copy
    my %primary_contig = Readonly::Clone %PRIMARY_CONTIG;

    ## No defined synonyms map
    return if ( not exists $primary_contig{$version}{synonyms_map} );

    $log->info( q{Creating VEP synonyms file: } . $outfile_path, $NEWLINE );

    ## Create dir if it does not exists
    make_path( dirname($outfile_path) );

    open my $FILEHANDLE_SYS, q{>}, $outfile_path
      or $log->logdie(qq{Cannot open $outfile_path: $ERRNO });

  CONTIG:
    foreach my $primary_contig ( @{ $primary_contig{$version}{contigs} } ) {

        ## Unpack
        my $synonymous_contig = $primary_contig{$version}{synonyms_map}{$primary_contig};

        ## No defined synonym
        next CONTIG if ( not $synonymous_contig );

        say {$FILEHANDLE_SYS} $primary_contig . $TAB . $synonymous_contig;
    }
    $log->info( q{Wrote: } . $outfile_path, $NEWLINE );
    close $FILEHANDLE_SYS;
    return $outfile_path;
}

1;
