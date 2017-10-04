package MIP::Log::MIP_log4perl;

use strict;
use warnings;
use warnings qw{ FATAL utf8 };
use utf8;
use open qw{ :encoding(UTF-8) :std };
use charnames qw{ :full :short };
use Carp;
use English qw{ -no_match_vars };
use Params::Check qw{ check allow last_error };

## CPANM
use Readonly;
use Log::Log4perl;

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.00;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw{ initiate_logger create_log4perl_congfig };
}

## Constants
Readonly my $COMMA => q{,};
Readonly my $NEWLINE => qq{\n};
Readonly my $SPACE => q{ };

sub initiate_logger {

## initiate_logger

## Function : Initiate the logger object
## Returns  : "$logger" {OBJ}
## Arguments: $categories_ref, $file_path, $log_name
##          : $categories_ref => Log categories {REF}
##          : $file_path      => log4perl config file path
##          : $log_name       => Log name

    my ($arg_href) = @_;

    ## Default(s)
    my $categories_ref;

    ## Flatten argument(s)
    my $file_path;
    my $log_name;

    my $tmpl = {
        categories_ref => {
            default     => [qw{ TRACE LogFile ScreenApp }],
            strict_type => 1,
            store       => \$categories_ref
        },
        file_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$file_path
        },
        log_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$log_name
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Creates log for the master script
    my $config = create_log4perl_congfig(
        {
            categories_ref => $categories_ref,
            file_path      => $file_path,
            log_name       => $log_name,
        }
    );

    Log::Log4perl->init( \$config );
    my $logger = Log::Log4perl->get_logger($log_name);
    return $logger;
}

sub create_log4perl_congfig {

## create_log4perl_congfig

## Function : Create log4perl config file.
## Returns  : $config
## Arguments: $categories_ref, $file_path, $log_name
##          : $categories_ref => Log categories {REF}
##          : $file_path      => log4perl config file path
##          : $log_name       => Log name

    my ($arg_href) = @_;

    ## Default(s)
    my $categories_ref;

    ## Flatten argument(s)
    my $file_path;
    my $log_name;

    my $tmpl = {
        categories_ref => {
            default     => [qw{ TRACE LogFile ScreenApp }],
            strict_type => 1,
            store       => \$categories_ref
        },
        file_path => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$file_path
        },
        log_name => {
            required    => 1,
            defined     => 1,
            strict_type => 1,
            store       => \$log_name
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $config = q{log4perl.category.} . $log_name . $SPACE . q{=} . $SPACE . join $COMMA . $SPACE, @{$categories_ref};
    $config .=<<"EOF";
$NEWLINE log4perl.appender.LogFile = Log::Log4perl::Appender::File
$NEWLINE log4perl.appender.LogFile.filename = $file_path
$NEWLINE log4perl.appender.LogFile.layout=PatternLayout
$NEWLINE log4perl.appender.LogFile.layout.ConversionPattern = [%p] %d %c - %m%n
$NEWLINE log4perl.appender.ScreenApp = Log::Log4perl::Appender::Screen
$NEWLINE log4perl.appender.ScreenApp.layout = PatternLayout
$NEWLINE log4perl.appender.ScreenApp.layout.ConversionPattern = [%p] %d %c - %m%n
EOF
    return $config;
}

1;
