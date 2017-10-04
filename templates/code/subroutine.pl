sub name_of_subroutine {

## name_of_subroutine

## Function :
## Returns  :

## Arguments: $hash_href, $arrays_ref, $scalar
##          : $hash_href  => Hash ref description {REF}
##          : $arrays_ref => Array ref description {REF}
##          : $scalar     => Scalar description

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $hash_href;
    my $arrays_ref;

    ## Default(s)
    my $scalar;

    my $tmpl = {
        hash_href => {
            required    => 1,
            defined     => 1,
            default     => {},
            strict_type => 1,
            store       => \$hash_href,
        },
        arrays_ref => {
            required    => 1,
            defined     => 1,
            default     => [],
            strict_type => 1,
            store       => \$arrays_ref,
        },
        scalar => {
            default     => 1,
            allow       => qr/ ^\d+$ /sxm,
            strict_type => 1,
            store       => \$scalar,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return;
}
