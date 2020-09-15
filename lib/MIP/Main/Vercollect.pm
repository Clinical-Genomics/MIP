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
use MIP::Constants qw{ $LOG_NAME $SPACE };
use MIP::Io::Read qw{ read_from_file };
use MIP::Io::Write qw{ write_to_file };

BEGIN {

    use base qw{ Exporter };
    require Exporter;

    # Set the version for version checking
    our $VERSION = q{1.0.4};

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ mip_vercollect };
}

sub mip_vercollect {

## Function : Execute mip vercollect to get executable versions
## Returns  :
## Arguments: $infile_path => Executables (binary) info file
##          : $outfile     => Data metric output file

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $outfile;
    my $infile_path;

    my $tmpl = {
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

    use MIP::Environment::Executable qw{ get_binaries_versions };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my %binary_info;

    ## Loads a YAML file into an arbitrary hash and returns it
    $log->info( q{Loading: } . $infile_path );
    my %binary_path = read_from_file(
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
        %binary_info = read_from_file(
            {
                format => q{yaml},
                path   => $outfile,
            }
        );
        $log->info( q{Loaded: } . $outfile );
    }

    ## Get/update executable versions
    %binary_info =
      ( %binary_info, get_binaries_versions( { binary_info_href => \%binary_path, } ) );

    ## Writes a binary hash to file
    write_to_file(
        {
            data_href => \%binary_info,
            format    => q{yaml},
            path      => $outfile,
        }
    );
    $log->info( q{Wrote: } . $outfile );

    return;
}

1;
