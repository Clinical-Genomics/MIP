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
use MIP::Constants qw{ $PIPE $SPACE };
use MIP::File::Format::Yaml qw{ load_yaml write_yaml };

BEGIN {

    use base qw{ Exporter };
    require Exporter;

    # Set the version for version checking
    our $VERSION = q{1.0.0};

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ get_binary version mip_vercollect };
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

    ## Loads a YAML file into an arbitrary hash and returns it
    $log->info( q{Loading: } . $infile_path );
    my %binary_info = load_yaml( { yaml_file => $infile_path, } );
    $log->info( q{Loaded: } . $infile_path );

    ## Get executable versions
    my %binary_version = get_binary_version( { binary_info_href => \%binary_info, } );

    ## Writes a qc data hash to file
    write_yaml(
        {
            yaml_file_path => $outfile . q{.yaml},
            yaml_href      => \%binary_version,
        }
    );
    $log->info( q{Wrote: } . $outfile . q{.yaml} );

    return;
}

######################
####Sub routines######
######################

sub get_binary_version {

## Function : Get executable versions
## Returns  : %binary_version
## Arguments: $binary_info_href         => Binary_Info_Href object

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $binary_info_href;

    my $tmpl = {
        binary_info_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$binary_info_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Get::Executable qw{ get_executable };

    my %binary_version;
    my %executable = get_executable( {} );

  BINARY:
    while ( my ( $binary, $binary_path ) = each %{$binary_info_href} ) {

        next BINARY if ( not exists $executable{$binary} );

        ## Unpack
        my $version_cmd    = $executable{$binary}{version_cmd};
        my $version_regexp = $executable{$binary}{version_regexp};

        my @get_version_cmds = _build_version_cmd(
            {
                binary_path    => $binary_path,
                version_cmd    => $version_cmd,
                version_regexp => $version_regexp,
            }
        );

        my %cmd_output = system_cmd_call(
            {
                command_string => join $SPACE,
                @get_version_cmds,
            }
        );
        ## Set binary version
        $binary_version{$binary} = $cmd_output{output}[0];
    }
    return %binary_version;
}

sub _build_version_cmd {

## Function : Execute mip vercollect to get executable versions
## Returns  : @cmds
## Arguments: $binary_path    => Executables (binary) file path
##          : $version_cmd    => Version command line option
##          : $version_regexp => Version reg exp to get version from system call

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $binary_path;
    my $version_cmd;
    my $version_regexp;

    my $tmpl = {
        binary_path => {
            defined     => 1,
            required    => 1,
            store       => \$binary_path,
            strict_type => 1,
        },
        version_cmd => {
            store       => \$version_cmd,
            strict_type => 1,
        },
        version_regexp => {
            defined     => 1,
            required    => 1,
            store       => \$version_regexp,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Language::Perl qw{ perl_nae_oneliners };
    use MIP::Unix::System qw{ system_cmd_call };

    ## Get perl wrapper around regexp
    my @perl_commands = perl_nae_oneliners(
        {
            oneliner_cmd => $version_regexp,
        }
    );

    my @get_version_cmds = ($binary_path);

    ## Allow for binary without any option to get version
    if ($version_cmd) {

        push @get_version_cmds, $version_cmd;
    }
    push @get_version_cmds, ( $PIPE, @perl_commands );
    return @get_version_cmds;
}

1;
