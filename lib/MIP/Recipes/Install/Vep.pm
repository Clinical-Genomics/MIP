package MIP::Recipes::Install::Vep;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Spec::Functions qw{ catdir catfile };
use List::MoreUtils qw{ any };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings qw{ FATAL utf8 };
use warnings;

## CPAN
use autodie qw{ :all };
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $BACKTICK $DASH $DOT $EQUALS $LOG_NAME $NEWLINE $PIPE $SPACE };
use MIP::Gnu::Coreutils qw{ gnu_mkdir gnu_rm };
use MIP::Language::Perl qw{ perl_nae_oneliners };
use MIP::Program::Tar qw{ tar };
use MIP::Program::Download::Wget qw{ wget };
use MIP::Program::Singularity qw{ singularity_exec };
use MIP::Program::Variantcalling::Vep qw{ variant_effect_predictor_install };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.20;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ install_vep };
}

sub install_vep {

## Function : Install plugins, download references and cache
## Returns  :
## Arguments: $active_parameter_href => Active parameter hash {REF}
##          : $container_href        => Container hash {REF}
##          : $container_path        => Path to VEP container
##          : $filehandle            => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $container_href;
    my $container_path;
    my $filehandle;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        container_href => {
            default     => {},
            required    => 1,
            store       => \$container_href,
            strict_type => 1,
        },
        container_href => {
            default     => {},
            required    => 1,
            store       => \$container_href,
            strict_type => 1,
        },
        container_path => {
            defined     => 1,
            required    => 1,
            store       => \$container_path,
            strict_type => 1,
        },
        filehandle => {
            defined  => 1,
            required => 1,
            store    => \$filehandle,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    my $auto_flag          = $active_parameter_href->{vep_auto_flag};
    my $cache_dir_path     = $active_parameter_href->{vep_cache_dir};
    my $reference_dir_path = $active_parameter_href->{reference_dir};
    my @assemblies         = @{ $active_parameter_href->{vep_assemblies} };
    my @plugins            = @{ $active_parameter_href->{vep_plugins} };
    my @species            = @{ $active_parameter_href->{vep_species} };

    ## Remove potential 'a' from auto flag since the API comes installed in the container
    $auto_flag =~ tr/a//d;

    ## Return if only API installation
    return if ( not $auto_flag );

    if ( not $cache_dir_path and not $reference_dir_path ) {
        $log->fatal(
q{Please supply a reference directory or a cache directory when installing VEP}
        );
        $log->fatal(
q{By default VEP cache and plugins will be downloaded to <reference_dir>/ensembl-tools-release-<version>/cache}
        );
        exit 1;
    }

    ## Install VEP
    say {$filehandle} q{## Install VEP plugins and cache};
    $log->info(q{Writing instructions for VEP installation});

    ## Get VEP version
    my $vep_version     = q?${VEP_VERSION}?;
    my $vep_version_cmd = _get_vep_version_cmd(
        {
            container_path => $container_path,
        }
    );
    say {$filehandle} q{VEP_VERSION} . $EQUALS . $vep_version_cmd;

    ## Setup cache_dir_path
    if ( not $cache_dir_path ) {
        $cache_dir_path = catdir( $reference_dir_path,
            q{ensembl-tools-release-} . $vep_version, q{cache} );
    }
    push @{ $container_href->{program_bind_paths} }, $cache_dir_path;

    ## Make sure that the cache directory exists
    if ( not -d $cache_dir_path ) {
        say {$filehandle} q{## Create cache directory};
        gnu_mkdir(
            {
                filehandle       => $filehandle,
                indirectory_path => $cache_dir_path,
                parents          => 1,
            }
        );
        say {$filehandle} $NEWLINE;
    }

    ## Don't install plugins unless specified in the auto flag
    if ( not $auto_flag =~ m/p/xms ) {
        undef @plugins;
    }

    my @vep_install_cmds = variant_effect_predictor_install(
        {
            assembly        => $assemblies[0],
            auto            => $auto_flag,
            cache_directory => $cache_dir_path,
            cache_version   => $vep_version,
            plugins_ref     => \@plugins,
            species_ref     => \@species,
        }
    );
    singularity_exec(
        {
            filehandle                     => $filehandle,
            singularity_container          => $container_path,
            singularity_container_cmds_ref => \@vep_install_cmds,
        }
    );
    say {$filehandle} $NEWLINE;

    ## If more than one assembly requested
    if ( ( scalar @assemblies > 1 ) && ( $auto_flag =~ / [cf] /xsm ) ) {

        ## Remove the plugins from the auto flag
        my $cf_auto_flag = $auto_flag;
        $cf_auto_flag =~ tr/p//d;

        # Find last index of array and initate
        Readonly my $NUMBER_OF_ASSEMBLIES => $#assemblies;

      ASSEMBLY:
        for my $assembly_version ( 1 .. $NUMBER_OF_ASSEMBLIES ) {
            ## Skip first assembly since it is already installed above

            say {$filehandle} q{## Install additional VEP cache assembly version};

            @vep_install_cmds = variant_effect_predictor_install(
                {
                    assembly        => $assemblies[$assembly_version],
                    auto            => $cf_auto_flag,
                    cache_directory => $cache_dir_path,
                    cache_version   => $vep_version,
                    species_ref     => \@species,
                }
            );
            singularity_exec(
                {
                    bind_paths_ref                 => [$cache_dir_path],
                    filehandle                     => $filehandle,
                    singularity_container          => $container_path,
                    singularity_container_cmds_ref => \@vep_install_cmds,
                }
            );
            say {$filehandle} $NEWLINE;
        }
    }

    ## Install and download extra plugin files
    if ( @plugins && $auto_flag =~ m/p/xms ) {

        my %finish_plugin_installation = (
            MaxEntScan => \&_install_maxentscan_plugin,
            LoFtool    => \&_install_loftool_plugin,
            ExACpLI    => \&_install_exacpli_plugin,
        );

      PLUGIN:
        foreach my $plugin (@plugins) {

            next PLUGIN if ( not $finish_plugin_installation{$plugin} );

            $finish_plugin_installation{$plugin}->(
                {
                    filehandle      => $filehandle,
                    plugin_dir_path => catdir( $cache_dir_path, q{Plugins} ),
                }
            );
        }
    }
    return 1;
}

sub _get_vep_version_cmd {

## Function : Write command that captures vep version in a bash variable
## Returns  : $vep_version_cmd
## Arguments: $container_path => Path to vep container

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $container_path;

    my $tmpl = {
        container_path => {
            defined     => 1,
            required    => 1,
            store       => \$container_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $vep_version_cmd = q{vep} . $SPACE . $PIPE . $SPACE;

    ## get perl oneliner to capture version number from output
    my @perl_commands = perl_nae_oneliners(
        {
            oneliner_name => q{get_vep_version},
        }
    );

    $vep_version_cmd .= join $SPACE, @perl_commands;

    my @vep_version_cmds = singularity_exec(
        {
            singularity_container          => $container_path,
            singularity_container_cmds_ref => [$vep_version_cmd],
        }
    );
    $vep_version_cmd = join $SPACE, @vep_version_cmds;

    return $BACKTICK . $vep_version_cmd . $BACKTICK;
}

sub _install_maxentscan_plugin {

## Function : Write command that installs MaxEntScan
## Returns  :
## Arguments: $filehandle      => filehandle
##          : $plugin_dir_path => Plugin path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $plugin_dir_path;
    my $filehandle;

    my $tmpl = {
        filehandle => {
            required => 1,
            store    => \$filehandle,
        },
        plugin_dir_path => {
            defined     => 1,
            required    => 1,
            store       => \$plugin_dir_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    say {$filehandle} q{## Add MaxEntScan required text file};
    my $maxent_file_path = catfile( $plugin_dir_path, q{fordownload.tar.gz} );
    wget(
        {
            filehandle   => $filehandle,
            outfile_path => $maxent_file_path,
            url =>
              q{http://hollywood.mit.edu/burgelab/maxent/download/fordownload.tar.gz},
        }
    );
    print {$filehandle} $NEWLINE;

    tar(
        {
            extract           => 1,
            outdirectory_path => dirname($maxent_file_path),
            filehandle        => $filehandle,
            filter_gzip       => 1,
            file_path         => $maxent_file_path,
        }
    );
    print {$filehandle} $NEWLINE;

    gnu_rm(
        {
            filehandle  => $filehandle,
            force       => 1,
            infile_path => $maxent_file_path,
        }
    );
    say {$filehandle} $NEWLINE;

    return;
}

sub _install_loftool_plugin {

## Function : Write command that downloads file required by LofTool
## Returns  :
## Arguments: $filehandle      => filehandle
##          : $plugin_dir_path => Plugin path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $plugin_dir_path;
    my $filehandle;

    my $tmpl = {
        filehandle => {
            required => 1,
            store    => \$filehandle,
        },
        plugin_dir_path => {
            defined     => 1,
            required    => 1,
            store       => \$plugin_dir_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    say {$filehandle} q{## Add LofTool required text file};
    wget(
        {
            filehandle   => $filehandle,
            outfile_path => catfile( $plugin_dir_path, q{LoFtool_scores.txt} ),
            url =>
q{https://raw.githubusercontent.com/Ensembl/VEP_plugins/master/LoFtool_scores.txt},
        }
    );
    say {$filehandle} $NEWLINE;

    return;
}

sub _install_exacpli_plugin {

## Function : Write command that downloads file required by ExACpLI
## Returns  :
## Arguments: $filehandle      => filehandle
##          : $plugin_dir_path => Plugin path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $plugin_dir_path;
    my $filehandle;

    my $tmpl = {
        filehandle => {
            required => 1,
            store    => \$filehandle,
        },
        plugin_dir_path => {
            defined     => 1,
            required    => 1,
            store       => \$plugin_dir_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    say {$filehandle} q{## Add pLI required value file};
    wget(
        {
            filehandle   => $filehandle,
            outfile_path => catfile( $plugin_dir_path, q{ExACpLI_values.txt} ),
            url =>
q{https://raw.githubusercontent.com/Ensembl/VEP_plugins/master/ExACpLI_values.txt},
        }
    );
    say {$filehandle} $NEWLINE;

    return;
}

1;
