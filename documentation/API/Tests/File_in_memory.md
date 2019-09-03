# Writing file to memory
Write a file to memory instead of disc:
```Perl
# For storing info to write
my $file_content;

## Store file content in memory by using referenced variable
open my $FILEHANDLE, q{>}, \$file_content
    or croak q{Cannot write to}
    . $SPACE
    . $file_content
    . $COLON
    . $SPACE
    . $OS_ERROR;

## Write to file
my $base_command = q{samtools};
say {$FILEHANDLE} $base_command;

## Close the filehandle
close $FILEHANDLE;
```

Then you can use this as variable to access the content of the file, for instance:
```Perl
my ($returned_base_command) = $file_content =~ /^($base_command)/ms;
```
