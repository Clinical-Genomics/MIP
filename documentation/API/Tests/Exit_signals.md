## Testing exit signals and log messages
Use Test::Trap when testing that sub routines produces an exit signal or warning.

```Perl  
## Load module
use Test::Trap;

### Testing exit
## Given a non-existing path
my $test_path = catfile( qw{ path to nowhere } );
## When trying to find the path
trap {
    path_finder( { test_path => $test_path, } )
};
## Then exit and throw FATAL log message
ok( $trap->exit, q{Exit if the path cannot be found} );
like( $trap->stderr, qr/FATAL/xms, q{Throw fatal log message} );

### Testing return and warn
## Given an ambiguous path
my $test_path = catfile( qw{ path to somewhere } );
## When trying to find the path
my @response = trap {
    path_finder( { test_path => $test_path, } )
};
## Then return the path and warn
is( $response[0], $test_path, q{Return ambiguous path} );
like( $trap->stderr, qr/WARN/xms, q{Throw warning} );

### Testing when croaking()
## Then exit and throw FATAL log message
is($trap->leaveby, q{die}, q{Exit if the path cannot be found});
like( $trap->die, qr/ERROR_MSG/xms, q{Throw error} );
```
