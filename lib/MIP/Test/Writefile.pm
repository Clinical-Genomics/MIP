package MIP::Test::Writefile;

use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };
use open qw{ :encoding(UTF-8) :std };
use strict;
use Test::More;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie;
use Readonly;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.02;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ test_write_to_file };
}

## Constants
Readonly my $COLON => q{:};
Readonly my $SPACE => q{ };

sub test_write_to_file {

##Function : Test of writing to file using anonymous FILEHANDLE
##Returns  :
##Arguments: $args_ref             => Arguments to function call
##         : $base_commands_ref    => Base commands {REF}
##         : $module_function_cref => Module method to test
##         : $separator            => Separator to use when writing

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $args_ref;
    my $base_commands_ref;
    my $module_function_cref;
    my $separator;

    my $tmpl = {
        args_ref => {
            default     => [],
            defined     => 1,
            required    => 1,
            store       => \$args_ref,
            strict_type => 1,
        },
        base_commands_ref => {
            defined     => 1,
			default => [],
            required    => 1,
            store       => \$base_commands_ref,
            strict_type => 1,
        },
        module_function_cref =>
          { defined => 1, required => 1, store => \$module_function_cref, },
        separator => {
            default     => q{ },
            store       => \$separator,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    # Create anonymous filehandle
    my $FILEHANDLE = IO::Handle->new();

    # Add new FILEHANDLE to args
    push @{$args_ref}, q{FILEHANDLE}, $FILEHANDLE;

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

    $module_function_cref->( { @{$args_ref} } );

    close $FILEHANDLE;

	my $base_command = join $separator, @{ $base_commands_ref };
    
	## Perform test
    my ($returned_base_command) = $file_content =~ /^($base_command)/ms;

    is( $returned_base_command, $base_command,
        q{Write commands with '} . $separator . q{' to file} );
    return;
}

1;
