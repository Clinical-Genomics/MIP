package MIP::Log::MIP_log4perl;

use 5.026;
use Carp;
use charnames qw{ :full :short };
use Cwd;
use English qw{ -no_match_vars };
use File::Spec::Functions qw{ catfile };
use File::Path qw{ make_path };
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ check allow last_error };
use strict;
use utf8;
use warnings;
use warnings qw{ FATAL utf8 };

## CPANM
use Log::Log4perl qw{ get_logger :levels };
use Readonly;

## MIPs lib/
use MIP::Constants
  qw{ $CLOSE_BRACKET $COMMA $DOT $EMPTY_STR $NEWLINE $OPEN_BRACKET $SPACE $TAB $UNDERSCORE };

BEGIN {
    require Exporter;
    use base qw{ Exporter };

    # Set the version for version checking
    our $VERSION = 1.07;

    # Functions and variables which can be optionally exported
    our @EXPORT_OK =
      qw{ create_log4perl_config get_default_log4perl_file get_log log_display_recipe_for_user initiate_logger retrieve_log };
}

sub create_log4perl_config {

## Function : Create log4perl config file
## Returns  : $config
## Arguments: $categories_ref => Log categories {REF}
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
            store       => \$categories_ref,
            strict_type => 1,
        },
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        log_name => {
            defined     => 1,
            required    => 1,
            store       => \$log_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $config =
      q{log4perl.category.} . $log_name . $SPACE . q{=} . $SPACE . join $COMMA . $SPACE,
      @{$categories_ref};
    $config .= <<"EOF";
$NEWLINE log4perl.appender.LogFile = Log::Log4perl::Appender::File
$NEWLINE log4perl.appender.LogFile.filename = $file_path
$NEWLINE log4perl.appender.LogFile.layout=PatternLayout
$NEWLINE log4perl.appender.LogFile.layout.ConversionPattern = [%p] %d %c - %m%n
$NEWLINE log4perl.appender.ScreenApp = Log::Log4perl::Appender::ScreenColoredLevels
$NEWLINE log4perl.appender.ScreenApp.layout = PatternLayout
$NEWLINE log4perl.appender.ScreenApp.layout.ConversionPattern = [%p] %d %c - %m%n
$NEWLINE log4perl.appender.ScreenApp.color.DEBUG=
$NEWLINE log4perl.appender.ScreenApp.color.INFO=
$NEWLINE log4perl.appender.ScreenApp.color.WARN=yellow
$NEWLINE log4perl.appender.ScreenApp.color.ERROR=red
$NEWLINE log4perl.appender.ScreenApp.color.FATAL=red
EOF
    return $config;
}

sub get_default_log4perl_file {

## Function : Set the default Log4perl file using supplied dynamic parameters
## Returns  : $log_file
## Arguments: $cmd_input       => User supplied info on cmd for log_file option {REF}
##          : $date            => The date
##          : $date_time_stamp => The date and time
##          : $outdata_dir     => Outdata directory
##          : $script          => The script that is executed

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $cmd_input;
    my $date;
    my $date_time_stamp;
    my $script;

    ## Default(s)
    my $outdata_dir;

    my $tmpl = {
        cmd_input => { store => \$cmd_input, strict_type => 1, },
        date      => {
            defined     => 1,
            required    => 1,
            store       => \$date,
            strict_type => 1,
        },
        date_time_stamp => {
            defined     => 1,
            required    => 1,
            store       => \$date_time_stamp,
            strict_type => 1,
        },
        script => {
            defined     => 1,
            required    => 1,
            store       => \$script,
            strict_type => 1,
        },
        outdata_dir => {
            default     => $arg_href->{outdata_dir} ||= getcwd(),
            store       => \$outdata_dir,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return $cmd_input if ( defined $cmd_input );

    ## No input from cmd i.e. create default logging directory and set default
    make_path( catfile( $outdata_dir, q{mip_log}, $date ) );

    ## Build log filename
    my $log_file = catfile( $outdata_dir, q{mip_log}, $date,
        $script . $UNDERSCORE . $date_time_stamp . $DOT . q{log} );

    ## Return default log file
    return $log_file;
}

sub get_log {

## Function : Create a log object, set log file in active_parameters and return log object
## Returns  : $log
## Arguments: $active_parameter_href => Active parameters for this analysis hash {REF}
##          : $date                  => Date
##          : $date_time_stamp       => Date and time stamp
##          : $log_name              => Log name
##          : $script                => The script that is executed

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $date;
    my $date_time_stamp;
    my $log_name;
    my $script;

    my $tmpl = {
        active_parameter_href => {
            default     => {},
            defined     => 1,
            required    => 1,
            store       => \$active_parameter_href,
            strict_type => 1,
        },
        date => {
            defined     => 1,
            required    => 1,
            store       => \$date,
            strict_type => 1,
        },
        date_time_stamp => {
            defined     => 1,
            required    => 1,
            store       => \$date_time_stamp,
            strict_type => 1,
        },
        log_name => {
            defined     => 1,
            required    => 1,
            store       => \$log_name,
            strict_type => 1,
        },
        script => {
            defined     => 1,
            required    => 1,
            store       => \$script,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    $active_parameter_href->{log_file} = get_default_log4perl_file(
        {
            cmd_input       => $active_parameter_href->{log_file},
            date            => $date,
            date_time_stamp => $date_time_stamp,
            outdata_dir     => $active_parameter_href->{outdata_dir},
            script          => $script,
        }
    );

    ## Creates log object
    my $log = initiate_logger(
        {
            file_path => $active_parameter_href->{log_file},
            log_name  => $log_name,
        }
    );
    return $log;
}

sub log_display_recipe_for_user {

## Function : Mutate recipe for display
## Returns  :
## Arguments: $indent_level => Number of TABS
##          : $log          => Log object
##          : $recipe       => Recipe

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $log;
    my $recipe;

    ## Defaults
    my $indent_level;

    my $tmpl = {
        indent_level => {
            default     => $EMPTY_STR,
            store       => \$indent_level,
            strict_type => 1,
        },
        log => {
            defined  => 1,
            required => 1,
            store    => \$log,
        },
        recipe => {
            defined     => 1,
            required    => 1,
            store       => \$recipe,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Adjust indent level
    if ($indent_level) {
        $indent_level = $TAB x $indent_level;
    }

    ## Replace "_" with $SPACE and make upper case
    my $display_recipe_name = uc join $SPACE, split /_/sxm, $recipe;
    $log->info( $indent_level . $OPEN_BRACKET . $display_recipe_name . $CLOSE_BRACKET );
    return;
}

sub initiate_logger {

## Function : Initiate the logger object
## Returns  : $logger {OBJ}
## Arguments: $categories_ref => Log categories {REF}
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
            store       => \$categories_ref,
            strict_type => 1,
        },
        file_path => {
            defined     => 1,
            required    => 1,
            store       => \$file_path,
            strict_type => 1,
        },
        log_name => {
            defined     => 1,
            required    => 1,
            store       => \$log_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    use MIP::Constants qw{ set_log_name_constant };

    ## Creates config for the log file
    my $config = create_log4perl_config(
        {
            categories_ref => $categories_ref,
            file_path      => $file_path,
            log_name       => $log_name,
        }
    );

    Log::Log4perl->init( \$config );
    my $logger = Log::Log4perl->get_logger($log_name);

    ## Set log name as constant
    set_log_name_constant( { log_name => $log_name, } );

    return $logger;
}

sub retrieve_log {

## Function  : Retrieves logger object and sets log level
## Returns   : $log
## Arguments : $level    => Set log level
##           : $log_name => Name of log
##           : $quiet    => Set log level to warn
##           : $verbose  => Set log level to debug

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $log_name;
    my $quiet;
    my $verbose;

    ## Default(s)
    my $level;

    my $tmpl = {
        level => {
            default => $INFO,
            allow   => [ $DEBUG, $INFO, $WARN, $ERROR, $WARN ],
            store   => \$level,
        },
        log_name => {
            defined     => 1,
            required    => 1,
            store       => \$log_name,
            strict_type => 1,
        },
        verbose => {
            allow => [ undef, 0, 1 ],
            store => \$verbose,
        },
        quiet => {
            allow => [ undef, 0, 1 ],
            store => \$quiet,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    ## Get logger
    my $log = get_logger($log_name);

    ## Set logger level
    if ($verbose) {
        $log->level($DEBUG);
    }
    elsif ($quiet) {
        $log->level($WARN);
    }
    else {
        $log->level($level);
    }
    return $log;
}

1;
