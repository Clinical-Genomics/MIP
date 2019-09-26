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
use MIP::File::Format::Yaml qw{ load_yaml write_yaml };

BEGIN {

    use base qw{ Exporter };
    require Exporter;

    # Set the version for version checking
    our $VERSION = q{1.0.0};

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

    use MIP::Get::Executable qw{ get_executable };

    ## Save final output data
    my %binary_version;

    ## Loads a YAML file into an arbitrary hash and returns it
    $log->info( q{Loading: } . $infile_path );
    my %binary_info = load_yaml( { yaml_file => $infile_path, } );
    $log->info( q{Loaded: } . $infile_path );

    my %executable = get_executable( {} );

    while ( my ( $binary, $binary_path ) = each %binary_info ) {

        say STDERR $binary . " " . $binary_path;

        if ( exists $executable{$binary} ) {

            say STDERR $executable{$binary}{version_regexp};
        }
    }

    ## Writes a qc data hash to file
    write_yaml(
        {
            yaml_file_path => $outfile,
            yaml_href      => \%binary_info,
        }
    );
    $log->info( q{Wrote: } . $outfile );

    return;
}

######################
####Sub routines######
######################

1;
