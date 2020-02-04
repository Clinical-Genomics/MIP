package MIP::Main::Vercollect;

#### Collects executable binary versions. Loads information on binaries to get versions
#### outputs exracted versions in YAML format.

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use File::Basename qw{ fileparse };
use File::Spec::Functions qw{ catdir catfile };
use FindBin qw{ $Bin };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ open close :all };
use Modern::Perl qw{ 2018 };

## MIPs lib/
use MIP::Constants qw{ $SPACE };
use MIP::Io::Read qw{ read_from_file write_to_file };

BEGIN {

    use base qw{ Exporter };
    require Exporter;

    # Set the version for version checking
    our $VERSION = q{1.0.1};

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ mip_vercollect };
}

sub mip_vercollect {

## Function : Execute mip vercollect to get executable versions
## Returns  :
## Arguments: $log         => Log object
##          : $outfile     => Data metric output file
##          : $infile_path => Executables (binary) info file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $log;
    my $outfile;
    my $infile_path;

    my $tmpl = {
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        outfile => {
            defined     => 1,
            required    => 1,
            store       => \$outfile,
            strict_type => 1,
        },
        infile_path => {
            defined     => 1,
            required    => 1,
            store       => \$infile_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::Executable qw{ get_binary_version };

    my %binary_version;

    ## Loads a YAML file into an arbitrary hash and returns it
    $log->info( q{Loading: } . $infile_path );
    my %binary_info = read_from_file(
        {
            format => q{yaml},
            path   => $infile_path,
        }
    );
    $log->info( q{Loaded: } . $infile_path );

    ## If outfile was previously generated - load to update
    if ( -e $outfile ) {
        ## Loads a YAML file into an arbitrary hash and returns it
        $log->info( q{Loading: } . $outfile );
        %binary_version = read_from_file(
            {
                format => q{yaml},
                path   => $outfile,
            }
        );
        $log->info( q{Loaded: } . $outfile );
    }

    ## Get/update executable versions
    %binary_version =
      ( %binary_version, get_binary_version( { binary_info_href => \%binary_info, } ) );

    ## Writes a qc data hash to file
    write_to_file(
        {
            data_href => \%binary_version,
            format    => q{yaml},
            path      => $outfile,
        }
    );
    $log->info( q{Wrote: } . $outfile );

    return;
}

1;
