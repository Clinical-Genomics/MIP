## Log in tests
When the sub routine you want to test requires a Log4perl object you can use this code recipe to quickly set up a test logger or just import an already prepared test log fixture using:

```
use MIP::Test::Fixtures qw{ test_log };
my $log = test_log( {} );
```

### Code recipe
We are going to need a temporary file to write the log config to. Use the core perl module File::Temp to initiate a temporary file that will be automatically deleted when the process that initiated the file exists.

```Perl
## Import the core module
use File::Temp;
```

Import, load and test the internal MIP module interface to Log4perl.

```Perl
# Import MIP::Log::MIP_log4perl and load initate_logger
use MIP::Log::MIP_log4perl qw{ initiate_logger };

## Modules with import
    my %perl_module = (
        q{MIP::Log::MIP_log4perl} => [qw{ initiate_logger }],
        q{MIP::Script::Utils}     => [qw{ help }],
    );

  PERL_MODULE:
    while ( my ( $module, $module_import ) = each %perl_module ) {
        use_ok( $module, @{$module_import} )
          or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
```

Create a test dir for writing log file

```Perl
## Create temp logger
my $test_dir = File::Temp->newdir();
my $test_log_path = catfile( $test_dir, q{test.log} );
```

Create the log object by naming the log and supplying the file to the sub initiate_logger.
```Perl
## Creates log object
my $log = initiate_logger(
    {
        file_path => $test_log_path,
        log_name  => q{TEST},
    }
);
```
