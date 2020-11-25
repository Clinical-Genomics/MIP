package MIP::Recipes::Install::Vep;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Basename qw{ dirname };
use File::Path qw{ make_path };
use File::Spec::Functions qw{ catdir catfile devnull };
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
use MIP::Constants qw{ $CONTAINER_MANAGER $LOG_NAME $PIPE $SPACE };
use MIP::Environment::Container qw{ run_container };
use MIP::Environment::Child_process qw{ child_process };
use MIP::Environment::Executable qw{ get_binary_version get_executable };
use MIP::Program::Tar qw{ tar };
use MIP::Program::Vep qw{ variant_effect_predictor_install };
use MIP::Program::Wget qw{ wget };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.26;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ install_vep };
}

sub install_vep {

## Function : Install plugins, download references and cache
## Returns  :
## Arguments: $active_parameter_href => Active parameter hash {REF}
##          : $container_href        => Container hash {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $container_href;

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
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Unpack parameters
    my $auto_flag          = $active_parameter_href->{vep_auto_flag};
    my $cache_dir_path     = $active_parameter_href->{vep_cache_dir};
    my $container_path     = $container_href->{uri};
    my $reference_dir_path = $active_parameter_href->{reference_dir};
    my @assemblies         = @{ $active_parameter_href->{vep_assemblies} };
    my @plugins            = @{ $active_parameter_href->{vep_plugins} };
    my @species            = @{ $active_parameter_href->{vep_species} };
    my $verbose            = $active_parameter_href->{verbose};

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
    $log->info(q{Installing VEP cache and plugins});

    ## Setup cache_dir_path
    $cache_dir_path = _setup_cache_dir(
        {
            active_parameter_href => $active_parameter_href,
            cache_dir_path        => $cache_dir_path,
            container_path        => $container_path,
            reference_dir_path    => $reference_dir_path,
            verbose               => $verbose,
        }
    );

    ## Don't install plugins unless specified in the auto flag
    if ( not $auto_flag =~ m/p/xms ) {

        undef @plugins;
    }

    my @vep_install_cmds = variant_effect_predictor_install(
        {
            assembly        => $assemblies[0],
            auto            => $auto_flag,
            cache_directory => $cache_dir_path,
            plugins_ref     => \@plugins,
            species_ref     => \@species,
        }
    );

    my @container_cmds = [
        run_container(
            {
                active_parameter_href => $active_parameter_href,
                bind_paths_ref        => [$cache_dir_path],
                container_path        => $container_path,
                container_manager     => $CONTAINER_MANAGER,
                container_cmds_ref    => \@vep_install_cmds,
            }
        )
    ];

    ## If more than one assembly requested
    if ( ( scalar @assemblies > 1 ) && ( $auto_flag =~ / [cf] /xsm ) ) {

        ## Remove the plugins from the auto flag
        my $cf_auto_flag = $auto_flag;
        $cf_auto_flag =~ tr/p//d;

        # Find last index of array and initiate
        Readonly my $NUMBER_OF_ASSEMBLIES => $#assemblies;

      ASSEMBLY:
        for my $assembly_version ( 1 .. $NUMBER_OF_ASSEMBLIES ) {
            ## Skip first assembly since it is already installed above

            @vep_install_cmds = variant_effect_predictor_install(
                {
                    assembly        => $assemblies[$assembly_version],
                    auto            => $cf_auto_flag,
                    cache_directory => $cache_dir_path,
                    species_ref     => \@species,
                }
            );

            push @container_cmds,
              [
                run_container(
                    {
                        active_parameter_href => $active_parameter_href,
                        bind_paths_ref        => [$cache_dir_path],
                        container_path        => $container_path,
                        container_manager  => $active_parameter_href->{container_manager},
                        container_cmds_ref => \@vep_install_cmds,
                    }
                )
              ];
        }
    }

    ## Install and download extra plugin files
    if ( @plugins && $auto_flag =~ m/p/xms ) {

        my $plugin_dir_path            = catdir( $cache_dir_path, q{Plugins} );
        my %finish_plugin_installation = (
            MaxEntScan => \&_install_maxentscan_plugin,
            LoFtool    => \&_install_loftool_plugin,
            ExACpLI    => \&_install_exacpli_plugin,
        );

      PLUGIN:
        foreach my $plugin (@plugins) {

            next PLUGIN if ( not $finish_plugin_installation{$plugin} );

            push @container_cmds,
              [
                $finish_plugin_installation{$plugin}->(
                    {
                        plugin_dir_path => $plugin_dir_path,
                    }
                )
              ];
        }
    }

  CONTAINER_CMDS_REF:
    foreach my $container_cmds_ref (@container_cmds) {

        my %process_return = child_process(
            {
                commands_ref => [ join $SPACE, @{$container_cmds_ref} ],
                process_type => q{ipc_cmd_run},
            }
        );
        if ( not $process_return{success} ) {

            $log->fatal(q{VEP installation failed});
            $log->logdie( $process_return{error_message} );
        }
    }
    return 1;
}

sub _install_maxentscan_plugin {

## Function : Build command that installs MaxEntScan
## Returns  : @cmds
## Arguments: $plugin_dir_path => Plugin path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $plugin_dir_path;

    my $tmpl = {
        plugin_dir_path => {
            defined     => 1,
            required    => 1,
            store       => \$plugin_dir_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @cmds = wget(
        {
            outfile_path => catfile( dirname( devnull() ), q{stdout} ),
            quiet        => 1,
            url =>
              q{http://hollywood.mit.edu/burgelab/maxent/download/fordownload.tar.gz},
            verbose => 0,
        }
    );
    push @cmds, $PIPE;

    push @cmds,
      tar(
        {
            extract           => 1,
            outdirectory_path => $plugin_dir_path,
            filter_gzip       => 1,
        }
      );

    return @cmds;
}

sub _install_loftool_plugin {

## Function : Build command that downloads file required by LofTool
## Returns  : @cmds
## Arguments: plugin_dir_path => Plugin path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $plugin_dir_path;

    my $tmpl = {
        plugin_dir_path => {
            defined     => 1,
            required    => 1,
            store       => \$plugin_dir_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @cmds = wget(
        {
            outfile_path => catfile( $plugin_dir_path, q{LoFtool_scores.txt} ),
            url =>
q{https://raw.githubusercontent.com/Ensembl/VEP_plugins/master/LoFtool_scores.txt},
        }
    );

    return @cmds;
}

sub _install_exacpli_plugin {

## Function : Write command that downloads file required by ExACpLI
## Returns  : @cmds
## Arguments: $plugin_dir_path => Plugin path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $plugin_dir_path;

    my $tmpl = {
        plugin_dir_path => {
            defined     => 1,
            required    => 1,
            store       => \$plugin_dir_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my @cmds = wget(
        {
            outfile_path => catfile( $plugin_dir_path, q{ExACpLI_values.txt} ),
            url =>
q{https://raw.githubusercontent.com/Ensembl/VEP_plugins/master/ExACpLI_values.txt},
        }
    );

    return @cmds;
}

sub _setup_cache_dir {

    ## Function : Setup vep cache directory
    ## Returns  : $cache_dir_path
    ## Arguments: $active_parameter_href => Active parameter hash {REF}
    ##          : $cache_dir_path     => Path to vep cache
    ##          : $container_path     => Path/uri to container
    ##          : $reference_dir_path => Path to reference directory
    ##          : $verbose            => Verbose

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $cache_dir_path;
    my $container_path;
    my $reference_dir_path;
    my $verbose;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        cache_dir_path => {
            required    => 1,
            store       => \$cache_dir_path,
            strict_type => 1,
        },
        container_path => {
            defined     => 1,
            required    => 1,
            store       => \$container_path,
            strict_type => 1,
        },
        reference_dir_path => {
            defined     => 1,
            required    => 1,
            store       => \$reference_dir_path,
            strict_type => 1,
        },
        verbose => {
            required    => 1,
            store       => \$verbose,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $log = Log::Log4perl->get_logger($LOG_NAME);

    if ( not $cache_dir_path ) {

        my %capture_version_cmd = get_executable(
            {
                executable_name => q{vep}
            }
        );

        my @vep_launch_cmds = run_container(
            {
                active_parameter_href => $active_parameter_href,
                container_cmds_ref    => [q{vep}],
                container_path        => $container_path,
                container_manager     => $CONTAINER_MANAGER,
            }
        );

        my $vep_version = get_binary_version(
            {
                binary                   => q{vep},
                binary_path              => join( $SPACE, @vep_launch_cmds ),
                capture_version_cmd_href => \%capture_version_cmd,
            }
        );
        if ( not $vep_version ) {

            $log->warn(q{Could not find vep version});
            $vep_version = q{unknown};
        }
        $cache_dir_path =
          catdir( $reference_dir_path, q{ensembl-tools-release-} . $vep_version,
            q{cache} );

        $log->warn( q{Setting VEP cache dir to: } . $cache_dir_path );
    }

    ## Make sure that the cache directory exists
    if ( not -d $cache_dir_path ) {

        make_path( $cache_dir_path, { verbose => $verbose } );
    }

    return $cache_dir_path;
}

1;
