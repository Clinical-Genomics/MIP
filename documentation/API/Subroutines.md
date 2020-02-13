# Sub routines
There are two main type of sub routines in MIP.
1. Sub routines that are imported/exported
2. Utility sub routines that only resides in the script or module they are created. These are named using an underscore in the beginning of the sub routine name i.e. `sub _a_utility_sub`

## Template
A template for sub routines are found in the [code dir].

The sub routine consists of four sections:

1. A documentation part with a mandatory header:
```Perl
## Function : Describe the sub routine function here
## Returns  : Name of variables returned. If none leave blank
## Arguments: $arrays_ref => Array ref description {REF}
##          : $hash_href  => Hash ref description {REF}
##          : $scalar     => Scalar description
```
2. Initilization of parameters:
```Perl
my ($arg_href) = @_;  #Always start with unpacking the parameters array

    ## Flatten argument(s)
    my $arrays_ref;  # Parameters without a default value
    my $hash_href;

    ## Default(s)  # This is for parameters with a default value
    my $scalar;
```

3. Checking the supplied parameters:
```Perl
my $tmpl = {
        arrays_ref => {
            default     => [], # Empty array ref
            defined     => 1, # Must be defined when passed
            required    => 1, # Must be supplied when passed
            store       => \$arrays_ref, # Where to store the parameter
            strict_type => 1, # Require correct type
        },
        hash_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$hash_href,
            strict_type => 1,
        },
        scalar => {
            allow       => qr/ ^\d+$ /sxm, # Set allowed value
            default     => 1,
            store       => \$scalar,
            strict_type => 1,
        },
    };
check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!}; # Check parameters according to template
```

4. The main part of the sub routine that actually does something.
5. Finally always end with a return statement
```Perl
return;
```
[code dir]: https://github.com/Clinical-Genomics/MIP/tree/master/templates/code/
