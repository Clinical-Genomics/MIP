# Writing file to memory
Write a file to memory instead of disc:
```Perl
# Create anonymous filehandle
my $FILEHANDLE = IO::Handle->new();

# For storing info to write
my $file_content;

## Store file content in memory by using referenced variable
open $FILEHANDLE, q{>}, \$file_content
    or croak q{Cannot write to}
    . $SPACE
    . $file_content
    . $COLON
    . $SPACE
    . $OS_ERROR;

# Close the filehandle
close $FILEHANDLE;
```

Then you can use this as variable to access the content of the file, for instance:
```Perl
my ($returned_base_command) = $file_content =~ /^($base_command)/ms;
```
