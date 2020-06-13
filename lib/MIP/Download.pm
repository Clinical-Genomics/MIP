package MIP::Download;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw{ :all };

## MIPs lib/
use MIP::Constants qw{ $COMMA };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ check_user_reference get_download_reference_attributes };
}

sub check_user_reference {

## Function : Check that the user supplied reference id and version
## Returns  : 1
## Arguments: $reference_genome_versions_ref => Reference genome build versions
##          : $reference_ref                 => Defined reference id and version
##          : $user_supplied_reference_ref   => User supplied reference id and version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $reference_genome_versions_ref;
    my $reference_ref;
    my $user_supplied_reference_ref;

    my $tmpl = {
        reference_genome_versions_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$reference_genome_versions_ref,
            strict_type => 1,
        },
        reference_ref => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$reference_ref,
            strict_type => 1,
        },
        user_supplied_reference_ref => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$user_supplied_reference_ref,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $log = Log::Log4perl->get_logger( uc q{mip_download} );

    my %has_no_file_attribute_for;

  REFERENCE:
    while ( my ( $reference_id, $versions_ref ) = each %{$user_supplied_reference_ref} ) {

        if ( not exists $reference_ref->{$reference_id} ) {

            $log->fatal( q{Cannot find reference key:} . $reference_id );
            exit 1;
        }

      REFERENCE_VERSION:
        foreach my $version ( @{$versions_ref} ) {

          GENOME_VERSION:
            foreach my $reference_genome_version ( @{$reference_genome_versions_ref} ) {

                my $has_attribute = get_download_reference_attributes(
                    {
                        genome_version => $reference_genome_version,
                        id             => $reference_id,
                        reference_href => $reference_ref,
                        version        => $version,
                    }
                );

                next REFERENCE if ($has_attribute);

                push @{ $has_no_file_attribute_for{$version} }, $reference_genome_version;
            }

            $log->warn(
                q{Cannot find version key: }
                  . $version
                  . q{ for reference key:}
                  . $reference_id
                  . q{ with genome build version: }
                  . join $COMMA,
                @{ $has_no_file_attribute_for{$version} },
            );
            $log->warn(qq{Skipping reference: $reference_id });
            next REFERENCE;
        }
    }
    return 1;
}

sub get_download_reference_attributes {

## Function : Get reference id version attributes per genome build for download
## Returns  : %{$reference_file_attribute_href}
## Arguments: $genome_version => Reference genome build versions
##          : $id             => Reference id
##          : $reference_href => Reference id and reference version per genome version
##          : $version        => Reference version

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $genome_version;
    my $id;
    my $reference_href;
    my $version;

    my $tmpl = {
        genome_version => {
            store       => \$genome_version,
            strict_type => 1,
        },
        id => {
            defined     => 1,
            required    => 1,
            store       => \$id,
            strict_type => 1,
        },
        reference_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$reference_href,
            strict_type => 1,
        },
        version => {
            store       => \$version,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use Data::Diver qw{ Dive };

    my $reference_file_attribute_href =
      Dive( $reference_href, ( $id, $genome_version, $version ) );

    return if ( not $reference_file_attribute_href );

    return %{$reference_file_attribute_href};
}

1;
