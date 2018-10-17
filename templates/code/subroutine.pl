sub name_of_subroutine {

## Function :
## Returns  :
## Arguments: $arrays_ref => Array ref description {REF}
##          : $hash_href  => Hash ref description {REF}
##          : $scalar     => Scalar description

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $arrays_ref;
    my $hash_href;

    ## Default(s)
    my $scalar;

    my $tmpl = {
        arrays_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$arrays_ref,
            strict_type => 1,
        },
        hash_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$hash_href,
            strict_type => 1,
        },
        scalar => {
            allow       => qr{ \A\d+\z }sxm,
            default     => 1,
            store       => \$scalar,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return;
}
