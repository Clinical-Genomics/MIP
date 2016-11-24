#!/usr/bin/env perl

###Will test perl modules and some selected funtions as well as vcf keys both in header and body. Adjusts dynamically according to supplied config file.
###Copyright 2016 Henrik Stranneheim

use Modern::Perl '2014';
use warnings qw( FATAL utf8 );
use autodie;
use v5.18;  #Require at least perl 5.18
use utf8;  #Allow unicode characters in this script
use open qw( :encoding(UTF-8) :std );
use charnames qw( :full :short );

use Params::Check qw[check allow last_error];
use File::Spec::Functions qw(catdir catfile devnull);

BEGIN {


    my @modules = ("YAML",
		   "Log::Log4perl",
		   "Path::Iterator::Rule",
		   "Set::IntervalTree",  # vcfParser
		   "Net::SSLeay",  # VEP
		   "LWP::Simple",  # VEP
		   "LWP::Protocol::https",  # VEP
		   "Archive::Zip",  # VEP
		   "DBI",  # VEP
		   "JSON",  # VEP
		   "DBD::mysql",  # VEP
		   "CGI",  # VEP
		   "Sereal::Encoder",  # VEP
		   "Sereal::Decoder",  # VEP
	);

    ## Evaluate that all modules required are installed
    eval_modules({modules_ref => \@modules,
		 });


    sub eval_modules {

	##eval_modules

	##Function : Evaluate that all modules required are installed
	##Returns  : ""
	##Arguments: $modules_ref
	##         : $modules_ref => Array of module names

	local $Params::Check::PRESERVE_CASE = 1;

	my ($arg_href) = @_;

	##Flatten argument(s)
	my $modules_ref;

	my $tmpl = {
	    modules_ref => { required => 1, default => [], strict_type => 1, store => \$modules_ref},
	};

	check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

	foreach my $module (@$modules_ref) {

	    $module =~s/::/\//g;  #Replace "::" with "/" since the automatic replacement magic only occurs for barewords.
	    $module .= ".pm";  #Add perl module ending for the same reason

	    eval {

		require $module;
	    };
	    if($@) {

		warn("NOTE: ".$module." not installed - Please install to run MIP.\n");
		warn("NOTE: Aborting!\n");
		exit 1;
	    }
	}
    }
}

use Test::More;
use Getopt::Long;
use FindBin qw($Bin); #Find directory of script
use File::Basename qw(dirname);
use File::Spec::Functions qw(catdir);

##Cpan
use YAML;

use vars qw($USAGE);

BEGIN {

    $USAGE =
	qq{test.t infile.vcf [VCF] config_file [YAML]
           -h/--help Display this help message
           -v/--version Display version
        };
}

my ($infile, $config_file);
my (%active_parameter, %vcfparser_data, %consequence_severity);

my $test_version = "2.0.0";

if(scalar(@ARGV) == 0) {

    print STDOUT $USAGE, "\n";
    exit;
}


############
####MAIN####
############

$infile = $ARGV[0];
$config_file = $ARGV[1];

###User Options
GetOptions('h|help' => sub { print STDOUT $USAGE, "\n"; exit;},  #Display help text
	   'v|version' => sub { print STDOUT "\ninstall.pl ".$test_version, "\n\n"; exit;},  #Display version number
    );

unless (defined($infile)) {

    print STDERR "Please supply an infile\n";
    exit;
}

unless (defined($config_file)) {

    print STDERR "Please supply a config file\n";
    exit;
}

## Test perl modules and functions
test_modules();

if (defined($config_file)) {  #Input from cmd

    ## Loads a YAML file into an arbitrary hash and returns it.
    %active_parameter = load_yaml($config_file);  #Load parameters from configfile
}


if($infile =~/.selected.vcf/) {

    ## Reads a file containg features to be annotated using range queries
    read_range_file({vcfparser_data_href => \%vcfparser_data,
		     range_coulumns_ref => \@{ $active_parameter{vcfparser_select_feature_annotation_columns} },
		     infile_path => catfile($active_parameter{reference_dir}, $active_parameter{vcfparser_select_file}),
		     range_file_key => "select_file",
		    });
}
else {  #Range file

    ## Reads a file containg features to be annotated using range queries
    read_range_file({vcfparser_data_href => \%vcfparser_data,
		     range_coulumns_ref => \@{ $active_parameter{vcfparser_range_feature_annotation_columns} },
		     infile_path => catfile($active_parameter{reference_dir}, $active_parameter{vcfparser_range_feature_file}),
		     range_file_key => "range_file",
		    });
}

## Reads infile in vcf format and parses annotations
read_infile_vcf({active_parameter_href => \%active_parameter,
		 vcfparser_data_href => \%vcfparser_data,
		 consequence_severity_href => \%consequence_severity,
		 infile => $infile,
		});

Test::More::done_testing();   # reached the end safely


######################
####SubRoutines#######
######################


sub test_modules {

##test_modules

##Function : Test perl modules and functions
##Returns  : ""
##Arguments:
##         :

    print STDOUT "\nTesting perl modules and selected functions\n\n";

    use FindBin qw($Bin); #Find directory of script

    ok(defined($Bin),"FindBin: Locate directory of script");

    use File::Basename qw(dirname);  #Strip the last part of directory

    ok(dirname($Bin), "File::Basename qw(dirname): Strip the last part of directory");

    use File::Spec::Functions qw(catdir);

    ok(catdir(dirname($Bin), "t"),"File::Spec::Functions qw(catdir): Concatenate directories");

    use YAML;

    my $yaml_file = catdir(dirname($Bin), "definitions", "define_parameters.yaml");
    ok( -f $yaml_file,"YAML: File= $yaml_file in MIP directory");

    my $yaml = YAML::LoadFile($yaml_file);  #Create an object
    ok( defined $yaml,"YAML: Load File" );  #Check that we got something
    ok(Dump( $yaml ),"YAML: Dump file");

    use Log::Log4perl;
    ## Creates log
    my $log_file = catdir(dirname($Bin), "templates", "mip_config.yaml");
    ok( -f $log_file,"Log::Log4perl: File= $log_file in MIP directory");

    ## Create log4perl config file
    my $config = create_log4perl_congfig(\$log_file);

    ok(Log::Log4perl->init(\$config), "Log::Log4perl: Initate");
    ok(Log::Log4perl->get_logger("mip_logger"), "Log::Log4perl: Get logger");

    my $logger = Log::Log4perl->get_logger("mip_logger");
    ok($logger->info("1"), "Log::Log4perl: info");
    ok($logger->warn("1"), "Log::Log4perl: warn");
    ok($logger->error("1"), "Log::Log4perl: error");
    ok($logger->fatal("1"), "Log::Log4perl: fatal");

    use Getopt::Long;
    push(@ARGV, ("-verbose", "2"));
    my $verbose = 1;
    ok(GetOptions("verbose:n"  => \$verbose), "Getopt::Long: Get options call");
    ok ($verbose == 2, "Getopt::Long: Get options modified");
}

sub create_log4perl_congfig {

##create_log4perl_congfig

##Function : Create log4perl config file.
##Returns  : "$config"
##Arguments: $file_name
##         : $file_name => log4perl config file {REF}

    my $file_name_ref = $_[0];

    my $conf = q?
        log4perl.category.mip_logger = TRACE, ScreenApp
        log4perl.appender.LogFile = Log::Log4perl::Appender::File
        log4perl.appender.LogFile.filename = ?.$$file_name_ref.q?
        log4perl.appender.LogFile.layout=PatternLayout
        log4perl.appender.LogFile.layout.ConversionPattern = [%p] %d %c - %m%n
        log4perl.appender.ScreenApp = Log::Log4perl::Appender::Screen
        log4perl.appender.ScreenApp.layout = PatternLayout
        log4perl.appender.ScreenApp.layout.ConversionPattern = [%p] %d %c - %m%n
        ?;
    return $conf;
}


sub read_infile_vcf {

##read_infile_vcf

##Function : Reads infile in vcf format and adds and parses annotations
##Returns  : ""
##Arguments: $active_parameter_href, $vcfparser_data_href, $consequence_severity_href, $infile
##         : $active_parameter_href     => The active parameters for this analysis hash {REF}
##         : vcfparser_data_href        => The keys from vcfParser i.e. range file and select file
##         : $consequence_severity_href => Consequence severity for SO-terms {REF}
##         : $infile                    => Infile

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $active_parameter_href;
    my $vcfparser_data_href;
    my $consequence_severity_href;
    my $infile;

    my $tmpl = {
	active_parameter_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$active_parameter_href},
	vcfparser_data_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$vcfparser_data_href},
	consequence_severity_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$consequence_severity_href},
	infile => { required => 1, defined => 1, strict_type => 1, store => \$infile},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my @vep_format_fields;
    my %vep_format_fields_column;

    my %meta_data;
    my %vcf_header;
    my %vcf_info_key;

    print STDOUT "\nTesting vcf header:\n\n";

    open(my $VCF, "<".$infile) or die "Cannot open ".$infile.":".$!, "\n";

    while (<$VCF>) {

	chomp $_;  # Remove newline
	if($. > 5000) {  #Exit after parsing X lines
	    last;
	}
	if (m/^\s+$/) {	# Avoid blank lines
	    next;
	}
	if ($_=~/^##(\S+)=/) {  # MetaData

	    parse_meta_data(\%meta_data, $_);

	    if ($_=~/INFO\=\<ID\=([^,]+)/) {

		$vcf_header{INFO}{$1} = $1; #Save to hash
	    }
	    if ($_=~/INFO\=\<ID\=CSQ/) { #Find VEP INFO Field

		ok($active_parameter_href->{pvarianteffectpredictor}, "VEP: CSQ in header and VEP should have been executed");

		if ($_=~/Format:\s(\S+)"\>/) { #Locate Format within VEP INFO meta line

		    @vep_format_fields = split(/\|/, $1);

		    for (my $field_counter=0;$field_counter<scalar(@vep_format_fields);$field_counter++) {

			$vep_format_fields_column{$vep_format_fields[$field_counter]} = $field_counter; #Save the order of VEP features
		    }
		}
		next;
	    }
	    next;
	}
	if ($_=~/^#CHROM/) {

	    ### Check Header now that we read all

	    ## VT
	    if ($active_parameter_href->{pvt} > 0) {

		if ($active_parameter_href->{vt_decompose} > 0) {

		    ok( defined($vcf_header{INFO}{OLD_MULTIALLELIC}), "VTDecompose key: OLD_MULTIALLELIC");
		}
		if ($active_parameter_href->{vt_normalize} > 0) {

		    ok( defined($vcf_header{INFO}{OLD_VARIANT}), "VTNormalize key: OLD_VARIANT");
		}
	    }

	    ## VCFParser
	    if ($active_parameter_href->{pvcfparser} > 0) {

		for my $key (keys %$vcfparser_data_href) {

		    ok( defined($vcf_header{INFO}{$key}), "VCFParser key: $key");
		}

		## Keys from vcfParser that are dynamically created from parsing the data
		my @vcfparser_dynamic_keys = ("most_severe_consequence",
		    );

		foreach my $key (@vcfparser_dynamic_keys) {

		    ok( defined($vcf_header{INFO}{$key}), "vcfparser dynamic keys: $key");
		}
	    }

	    ## Snpeff
	    if ($active_parameter_href->{psnpeff} > 0) {

		my @splitted_keys;
		for my $annotation_file (keys %{ $active_parameter_href->{snpsift_annotation_files} }) {

		    if ($active_parameter_href->{snpsift_annotation_outinfo_key}{$annotation_file}) {

			@splitted_keys = split(',', $active_parameter_href->{snpsift_annotation_outinfo_key}{$annotation_file});

			my @original_vcf_keys = split(',', $active_parameter_href->{snpsift_annotation_files}{$annotation_file});

			## Modify list elements in place to produce -names flag from SnpEff
			map{$_ .= shift(@original_vcf_keys)} @splitted_keys;
		    }
		    else {

			@splitted_keys = split(',', $active_parameter_href->{snpsift_annotation_files}{$annotation_file});
		    }
		    foreach my $key (@splitted_keys) {

			ok( defined($vcf_header{INFO}{$key}), "Snpsift annotation key: $key");
		    }
		}
		for my $key (@{ $active_parameter_href->{snpsift_dbnsfp_annotations} }) {

		    $key =~ s/\+/_/g;  #Special case due to the fact that snpEff v4.2 transforms + to _ for some reason
		    $key = "dbNSFP_".$key;
		    ok( defined($vcf_header{INFO}{$key}), "Snpsift dbNSFP_key: $key");
		}
	    }
	    ## GENMOD
	    ## 1000G key from genmodFilter
	    if ( ($active_parameter_href->{pvt} > 0) && ($active_parameter_href->{vt_genmod_filter} > 0) ) {

		if (defined($active_parameter_href->{vt_genmod_filter_1000G})) {

		    ok( defined($vcf_header{INFO}{'1000GAF'}), "Genmod filter: 1000GAF key");
		}
	    }
	    if ($active_parameter_href->{prankvariant} > 0) {

		## Keys from genmod
		my @genmod_keys = ("Compounds",
				   "RankScore",
				   "ModelScore",
				   "GeneticModels",
		    );

		foreach my $key (@genmod_keys) {

		    ok( defined($vcf_header{INFO}{$key}), "Genmod: $key");
		}

		## Spidex key from genmodAnnotate
		if ($active_parameter_href->{spidex_file}) {

		    ok( defined($vcf_header{INFO}{SPIDEX}), "Genmod annotate: SPIDEX key");
		}
		## CADD key from genmodAnnotate
		if ( ($active_parameter_href->{genmod_annotate_cadd_files}) || ($active_parameter_href->{genmod_annotate_cadd_files}) ) {

		    ok( defined($vcf_header{INFO}{CADD}), "GENMODAnnotate: CADD key");
		}
	    }
	    next;
	}
	if ( $_ =~/^(\S+)/ ) {

	    my @line_elements = split("\t",$_); #Loads vcf elements

	    my @key_values = split(/;/, $line_elements[7]); #Split INFO field to key=value items

	    for my $element (@key_values) {

		my @keys = split("=", $element);  #key = 0 and value = 1

		$vcf_info_key{$keys[0]}++;  #Increment
	    }
	}
    }
    close($VCF);

    ##Check keys found in INFO field

    print STDOUT "\nTesting vcf INFO fields and presence in header:\n\n";

    foreach my $key (keys %vcf_info_key) {

	ok( defined($vcf_header{INFO}{$key}), "Found both header and line field key for: $key. key count: ".$vcf_info_key{$key});
    }
}


sub read_range_file {

##read_range_file

##Function : Reads a file containg features to be annotated using range queries e.g. EnsemblGeneID.
##Returns  : ""
##Arguments: $vcfparser_data_href, $range_coulumns_ref, $range_file_key, $infile_path, $paddingRef, $select_feature_matching_column
##         : $vcfparser_data_href            => Range file hash {REF}
##         : $range_coulumns_ref             => Range columns to include {REF}
##         : $range_file_key                 => Range file key used to seperate range file(s) i.e., select and range
##         : $infile_path                    => Infile path
##         : $select_feature_matching_column => Column in the select file to match with vcf key annotation

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $vcfparser_data_href;
    my $range_coulumns_ref;
    my $range_file_key;
    my $infile_path;
    my $select_feature_matching_column;

    my $tmpl = {
	vcfparser_data_href => { required => 1, defined => 1, default => {}, strict_type => 1, store => \$vcfparser_data_href},
	range_coulumns_ref => { required => 1, defined => 1, default => [], strict_type => 1, store => \$range_coulumns_ref},
	range_file_key => { required => 1, defined => 1, strict_type => 1, store => \$range_file_key},
	infile_path => { required => 1, defined => 1, strict_type => 1, store => \$infile_path},
	select_feature_matching_column => { strict_type => 1, store => \$select_feature_matching_column},
    };

    check($tmpl, $arg_href, 1) or die qw[Could not parse arguments!];

    my @headers; #Save headers from rangeFile

    open(my $RRF, "<".$infile_path) or die "Cannot open ".$infile_path.":".$!, "\n";

    while (<$RRF>) {

	chomp $_; #Remove newline

	if (m/^\s+$/) {		# Avoid blank lines

	    next;
	}
	if ($_=~/^##/) {#MetaData - Avoid

	    next;
	}
	if ($_=~/^#/) {#Header/Comment

	    @headers = split(/\t/, $_);

	    for (my $extract_columns_counter=0;$extract_columns_counter<scalar(@$range_coulumns_ref);$extract_columns_counter++) { #Defines what scalar to store

		my $header_key_ref = \$headers[ $$range_coulumns_ref[$extract_columns_counter] ];
		$vcfparser_data_href->{$$header_key_ref} = $extract_columns_counter; #Column position in supplied range input file
	    }
	    next;
	}
    }
    close($RRF);
    print STDOUT "\nFinished reading ".$range_file_key." file: ".$infile_path,"\n";
}


sub parse_meta_data {

##parse_meta_data

##Function : Writes metadata to filehandle specified by order in meta_data_orders.
##Returns  : ""
##Arguments: $meta_data_href, $meta_data_string
##         : $meta_data_href   => Hash for meta_data {REF}
##         : $meta_data_string => The meta_data string from vcf header

    my $meta_data_href = $_[0];
    my $meta_data_string = $_[1];

    if ($meta_data_string=~/^##fileformat/) {  #Catch fileformat as it has to be at the top of header

	push(@{ $meta_data_href->{fileformat}{fileformat} }, $meta_data_string);  #Save metadata string
    }
    elsif ($meta_data_string=~/^##contig/) {  #catch contigs to not sort them later

	push(@{ $meta_data_href->{contig}{contig} }, $meta_data_string);  #Save metadata string
    }
    elsif ($meta_data_string=~/^##(\w+)=(\S+)/) {  #FILTER, FORMAT, INFO etc and more custom records

	push(@{ $meta_data_href->{$1}{$2} }, $meta_data_string);  #Save metadata string
    }
    else {  #All oddities

	push(@{ $meta_data_href->{other}{other} }, $meta_data_string);  #Save metadata string
    }
}


sub load_yaml {

##load_yaml

##Function : Loads a YAML file into an arbitrary hash and returns it. Note: Currently only supports hashreferences and hashes and no mixed entries.
##Returns  : %yaml_hash
##Arguments: $yaml_file
##         : $yaml_file => The yaml file to load

    my $yaml_file = $_[0];

    my %yaml_hash;

    open (my $YAML, "<", $yaml_file) or die "cannot open ".$yaml_file.":".$!, "\n";  #Log4perl not initialised yet, hence no logdie

    local $YAML::QuoteNumericStrings=1;  #Force numeric values to strings in YAML representation
    if (defined ($YAML::QuoteNumericStrings) && $YAML::QuoteNumericStrings == 1) {

	%yaml_hash = %{ YAML::LoadFile($yaml_file) };  #Load hashreference as hash
    }
    close($YAML);

    return %yaml_hash;
}
