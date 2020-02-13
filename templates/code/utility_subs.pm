## Utility subroutines that can be added to modules

sub _this_sub {

## Function : Returns the name of the current subroutine
## Returns  : $this_sub
## Arguments:

    ## Get full path to subroutine
    my $this_sub = ( caller 1 )[$THREE];

    ## Isolate subroutine
    $this_sub = ( split /::/xms, $this_sub )[$MINUS_ONE];

    return $this_sub;
}

sub _parent_module {

## Function : Returns the name of the module that called this one
## Returns  : $parent_module
## Arguments:

    ## Get full path to module
    my $parent_module = ( caller 1 )[0];

    ## Isolate module
    $parent_module = ( split /::/xms, $parent_module )[$MINUS_ONE];

    return $parent_module;
}

sub _parent_sub {

## Function : Returns the name of the parent subroutine
## Returns  : $parent_sub
## Arguments:

    ## Get full path to parent subroutine
    my $parent_sub = ( caller 2 )[$THREE];

    ## Isolate subroutine
    $parent_sub = ( split /::/xms, $parent_sub )[$MINUS_ONE];

    return $parent_sub;
}

