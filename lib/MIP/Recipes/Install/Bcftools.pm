package MIP::Recipes::Install::Bcftools;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd q{abs_path};
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
use MIP::Constants qw{ $DOUBLE_QUOTE };
use MIP::Gnu::Coreutils qw{ gnu_mkdir };
use MIP::Program::Singularity qw{ singularity_exec };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ install_bcftools };
}

sub install_bcftools {

## Function : Add relative bind path for bcftools
## Returns  :
## Arguments: $active_parameter_href => Active parameter hash {REF}
##          : $contaienr_href        => Container hah {REF}
##          : $container_path        => Path to container
##          : $filehandle            => Filehandle to write to

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $container_path;
    my $container_href;
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
        container_path => {
            defined     => 1,
            required    => 1,
            store       => \$container_path,
            strict_type => 1,
        },
        filehandle => {
            defined => 1,
            store   => \$filehandle,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Only add MIP_BIND for bcftools
    return if ( not $container_path =~ /bcftools.sif \z/sxm );

    my $bcftools_bind_path = $DOUBLE_QUOTE . q{$MIP_BIND} . $DOUBLE_QUOTE;

    ## Store annotation dir path for later
    if ( $container_href->{program_bind_paths} ) {

        push @{ $container_href->{program_bind_paths} }, $bcftools_bind_path;
    }
    else {
        $container_href->{program_bind_paths} = [$bcftools_bind_path];
    }
    return 1;
}

1;
