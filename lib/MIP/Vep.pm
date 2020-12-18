package MIP::Vep;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catdir catfile };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };
use List::MoreUtils qw { any };

## MIPs lib/
use MIP::Constants qw{ %CONTAINER_CMD $LOG_NAME $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{
      check_vep_api_cache_versions
      check_vep_custom_annotation
      check_vep_plugin
    };
}

sub check_vep_api_cache_versions {

## Function : Compare VEP API and VEP chache versions. Exit if non-match
## Returns  :
## Arguments: $vep_directory_cache => VEP cache directory path {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $vep_directory_cache;

    my $tmpl = {
        vep_directory_cache => {
            defined     => 1,
            required    => 1,
            store       => \$vep_directory_cache,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use File::Find::Rule;
    use MIP::Environment::Executable qw{ get_binary_version get_executable };

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    my %capture_version_cmd = get_executable( { executable_name => q{vep} } );
    my $vep_version         = get_binary_version(
        {
            binary                   => q{vep},
            binary_path              => $CONTAINER_CMD{vep},
            capture_version_cmd_href => \%capture_version_cmd,
        }
    );

    ## Check that a version number was picked up
    if ( not $vep_version ) {
        $log->warn(
            q{Could not retrieve VEP version. Skipping checking that VEP api and cache matches.});
        return;
    }

    ## Get absolute path to VEP species cache as it is commonly linked
    my $vep_cache_species_dir_path =
      _get_vep_cache_species_dir_path( { vep_directory_cache => $vep_directory_cache, } );

    ## Get folders in cache directory
    # Build rule
    my $rule = File::Find::Rule->new();

    # Find directories
    $rule->directory;

    # Only get directories in the current folder
    $rule->maxdepth(1);

    # Get relative paths
    $rule->relative;

    # Apply rule to find directories
    my @vep_cache_version_folders = $rule->in($vep_cache_species_dir_path);

    ## Check that
    if ( not @vep_cache_version_folders ) {
        $log->warn(
q{Could not retrieve VEP cache version. Skipping checking that VEP api and cache matches.}
        );
        return;
    }

    ## Check if the VEP api version and cache versions matches
    if ( not any { /$vep_version/xms } @vep_cache_version_folders ) {

        $log->fatal( q{Differing versions between 'VEP API':}
              . $SPACE
              . $vep_version
              . $SPACE
              . q{and '--vep_directory_cache':}
              . $SPACE
              . $vep_directory_cache );
        exit 1;
    }
    return 1;
}

sub check_vep_custom_annotation {

## Function : Check VEP custom annotations options
## Returns  :
## Arguments: $vep_custom_ann_href => VEP custom annotation {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $vep_custom_ann_href;

    my $tmpl = {
        vep_custom_ann_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$vep_custom_ann_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Path qw { check_filesystem_objects_and_index_existance };

    ## Nothing to check
    return 0 if ( not keys %{$vep_custom_ann_href} );

  ANN:
    while ( my ( $ann, $value_href ) = each %{$vep_custom_ann_href} ) {

        my $err_msg = $ann . q{ Is not a hash ref for vep_custom_annotation};
        croak($err_msg) if ( ref $value_href ne q{HASH} );

        ## Check the VEP custom annotations options and that they have allowed values
        _check_vep_custom_annotation_options(
            {
                annotation             => $ann,
                custom_ann_option_href => $value_href,
            }
        );

        ## Check path object exists
        check_filesystem_objects_and_index_existance(
            {
                object_name    => $ann,
                object_type    => q{file},
                parameter_name => q{vep_custom_annotation},
                path           => $value_href->{path},
            }
        );
    }
    return 1;
}

sub check_vep_plugin {

## Function : Check VEP plugins and options
## Returns  : 0 or 1
## Arguments: $parameter_name       => Parameter name
##          : $vep_plugin_href      => VEP plugin annotation {REF}
##          : $vep_plugins_dir_path => VEP plugin dir path

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $parameter_name;
    my $vep_plugin_href;
    my $vep_plugins_dir_path;

    my $tmpl = {
        parameter_name => {
            defined  => 1,
            required => 1,
            store    => \$parameter_name,
        },
        vep_plugin_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$vep_plugin_href,
            strict_type => 1,
        },
        vep_plugins_dir_path => {
            defined     => 1,
            required    => 1,
            store       => \$vep_plugins_dir_path,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::File::Path qw { check_filesystem_objects_and_index_existance };

    ## Nothing to check
    return 0 if ( not keys %{$vep_plugin_href} );

  PLUGIN:
    while ( my ( $plugin, $value_href ) = each %{$vep_plugin_href} ) {

        my $err_msg = $plugin . q{ Is not a hash ref for vep_plugin};
        croak($err_msg) if ( ref $value_href ne q{HASH} );

        check_filesystem_objects_and_index_existance(
            {
                object_name    => $plugin,
                object_type    => q{file},
                parameter_name => $parameter_name,
                path           => catfile( $vep_plugins_dir_path, $plugin . q{.pm} ),
            }
        );

        next PLUGIN if ( not exists $value_href->{exist_check} );

      OBJECT:
        foreach my $object_href ( @{ $value_href->{exist_check} } ) {

            check_filesystem_objects_and_index_existance(
                {
                    object_name    => $plugin,
                    object_type    => $object_href->{type},
                    parameter_name => $parameter_name,
                    path           => $object_href->{path},
                }
            );
        }
    }
    return 1;
}

sub _get_vep_cache_species_dir_path {

## Function : Get the vep cache species dir path
## Returns  : $vep_cache_species_dir_path
## Arguments: $vep_directory_cache => VEP cache directory path {REF}

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $vep_directory_cache;

    my $tmpl = {
        vep_directory_cache => {
            defined     => 1,
            required    => 1,
            store       => \$vep_directory_cache,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Cwd qw{ abs_path };

    my @vep_species_cache = qw{ homo_sapiens homo_sapiens_merged };
    my $vep_cache_species_dir_path;

  SPECIES_CACHE:
    foreach my $species_cache (@vep_species_cache) {

        ## Get species specific cache dir
        $vep_cache_species_dir_path = abs_path( catdir( $vep_directory_cache, $species_cache ) );
        last if ( -e $vep_cache_species_dir_path );
    }
    return $vep_cache_species_dir_path;
}

sub _check_vep_custom_annotation_options {

## Function : Check the VEP custom annotations options are defined and with allowed values
## Returns  :
## Arguments: $annotation             => Annotation
##          : $custom_ann_option_href => Custom annotation options

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $annotation;
    my $custom_ann_option_href;

    my $tmpl = {
        annotation => {
            defined     => 1,
            required    => 1,
            store       => \$annotation,
            strict_type => 1,
        },
        custom_ann_option_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$custom_ann_option_href,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger($LOG_NAME);

    ## Check required keys
    my @required_options = (qw{ key });
  REQ_OPTION:
    foreach my $required_option (@required_options) {

        if (    not exists $custom_ann_option_href->{$required_option}
            and not defined $custom_ann_option_href->{$required_option} )
        {

            $log->fatal( q{Vep custom annotation option hash: }
                  . $annotation
                  . q{ lacks required option }
                  . $required_option );
            exit 1;
        }
    }

    ## Check allowed and defined options for annotation
    my %check_vep_annotations = (
        annotation_type          => { allow   => [qw{ exact overlap }], },
        file_type                => { allow   => [qw{ bed gff gtf vcf bigwig }], },
        force_report_coordinates => { allow   => [ 0, 1 ], },
        key                      => { defined => 1, },
        path                     => { defined => 1, },
        vcf_fields               => { defined => 1, },
    );

  OPTION:
    foreach my $option ( keys %{$custom_ann_option_href} ) {

        ## Allow anything defined
        next OPTION if ( $check_vep_annotations{$option}{defined} );

        next OPTION
          if (
            any { $_ eq $custom_ann_option_href->{$option} }
            @{ $check_vep_annotations{$option}{allow} }
          );

        $log->fatal( q{Vep custom annotation option hash: }
              . $annotation
              . q{ has a not allowed option value '}
              . $option . q{ => }
              . $custom_ann_option_href->{$option} );
        $log->fatal( q{Allowed options are: } . join $SPACE,
            @{ $check_vep_annotations{$option}{allow} } );
        exit 1;
    }
    return 1;
}

1;
