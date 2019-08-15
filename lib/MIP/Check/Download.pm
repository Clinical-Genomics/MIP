package MIP::Check::Download;

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
use Readonly;

## MIPs lib/
use MIP::Constants qw{ $COMMA $SPACE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ check_user_reference };
}

sub check_user_reference {

## Function : Check that the user supplied reference id and version
## Returns  :
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

    ## Retrieve logger object
    my $log = Log::Log4perl->get_logger( uc q{mip_download} );

    ## Store what has been seen
    my %cache;

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

                ## Found match
                if (
                    exists $reference_ref->{$reference_id}{$reference_genome_version}
                    {$version} )
                {

                    $cache{$reference_id}++;
                }
                ## Store mismatch
                push @{ $cache{$version} }, $reference_genome_version;
            }

            ## Require at least one match
            next REFERENCE if ( $cache{$reference_id} );

            $log->warn(
                q{Cannot find version key: }
                  . $version
                  . q{ for reference key:}
                  . $reference_id
                  . q{ with genome build version: }
                  . join $COMMA,
                @{ $cache{$version} },
            );
$log->warn(qq{Skipping reference: $reference_id });
            next REFERENCE;
        }

    }
    return 1;
}

1;
