package MIP_log::Log4perl;

use Modern::Perl '2014';
use warnings qw( FATAL utf8 );
use autodie;
use v5.18;  #Require at least perl 5.18
use utf8;  #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );

BEGIN {
    require Exporter;

    # Set the version for version checking
    our $VERSION = 1.00;

    # Inherit from Exporter to export functions and variables
    our @ISA = qw(Exporter);

    # Functions and variables which are exported by default
    our @EXPORT = qw();

    # Functions and variables which can be optionally exported
    our @EXPORT_OK = qw(initiate_logger create_log4perl_congfig);

}

use Params::Check qw[check allow last_error];
$Params::Check::PRESERVE_CASE = 1;  #Do not convert to lower case

## Third party module(s)
use Log::Log4perl;


sub initiate_logger {
    
##initiate_logger
    
##Function : Initiate the logger object
##Returns  : "$logger" {OBJ}
##Arguments: $categories_ref, $file_path_ref, $log_name
##         : $categories_ref => Log categories {REF}
##         : $file_path_ref  => log4perl config file path {REF}
##         : $log_name       => Log name
    
    my ($arg_href) = @_;
    
    ## Default(s)
    my $categories_ref;

    ## Flatten argument(s)
    my $file_path_ref;
    my $log_name;

    my $tmpl = {
	categories_ref => { default => ["TRACE", "LogFile", "ScreenApp"],
			    strict_type => 1, store => \$categories_ref},
	file_path_ref => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$file_path_ref},
	log_name => { required => 1, defined => 1, strict_type => 1, store => \$log_name},
    };
    
    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];
    
    ## Creates log for the master script
    my $config = create_log4perl_congfig({categories_ref => $categories_ref,
					  file_path_ref => $file_path_ref,
					  log_name => $log_name,
					 });
    
    Log::Log4perl->init(\$config);
    my $logger = Log::Log4perl->get_logger($log_name);
    return $logger;
}


sub create_log4perl_congfig {

##create_log4perl_congfig

##Function : Create log4perl config file.
##Returns  : "$config"
##Arguments: $categories_ref, $file_path_ref, $log_name
##         : $categories_ref => Log categories {REF}
##         : $file_path_ref  => log4perl config file path {REF}
##         : $log_name       => Log name

    my ($arg_href) = @_;

    ## Default(s)
    my $categories_ref;

    ## Flatten argument(s)
    my $file_path_ref;
    my $log_name;

    my $tmpl = {
	categories_ref => { default => ["TRACE", "LogFile", "ScreenApp"],
			    strict_type => 1, store => \$categories_ref},
	file_path_ref => { required => 1, defined => 1, default => \$$, strict_type => 1, store => \$file_path_ref},
	log_name => { required => 1, defined => 1, strict_type => 1, store => \$log_name},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];
    my $config = q?
        log4perl.category.?.$log_name.q? = ?.join(", ", @$categories_ref).q?
        log4perl.appender.LogFile = Log::Log4perl::Appender::File
        log4perl.appender.LogFile.filename = ?.$$file_path_ref.q?
        log4perl.appender.LogFile.layout=PatternLayout
        log4perl.appender.LogFile.layout.ConversionPattern = [%p] %d %c - %m%n

        log4perl.appender.ScreenApp = Log::Log4perl::Appender::Screen
        log4perl.appender.ScreenApp.layout = PatternLayout
        log4perl.appender.ScreenApp.layout.ConversionPattern = [%p] %d %c - %m%n
        ?;
    return $config;
}


1;
