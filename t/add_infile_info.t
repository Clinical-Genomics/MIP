#!/usr/bin/env perl

use 5.018;
use Carp;
use charnames qw{ :full :short };
use English qw{ -no_match_vars };
use File::Basename qw{ basename dirname };
use File::Spec::Functions qw{ catdir catfile };
use File::Temp;
use FindBin qw{ $Bin };
use Getopt::Long;
use open qw{ :encoding(UTF-8) :std };
use Params::Check qw{ allow check last_error };
use Test::More;
use Time::Piece;
use utf8;
use warnings qw{ FATAL utf8 };

## CPANM
use autodie qw { :all };
use Modern::Perl qw{ 2014 };
use Readonly;

## MIPs lib/
use lib catdir( dirname($Bin), q{lib} );
use MIP::Log::MIP_log4perl qw{ initiate_logger };
use MIP::Script::Utils qw{ help };

our $USAGE = build_usage( {} );

my $VERBOSE = 1;
our $VERSION = '1.0.0';

## Constants
Readonly my $COMMA       => q{,};
Readonly my $DOT         => q{.};
Readonly my $NEWLINE     => qq{\n};
Readonly my $READ_LENGTH => 151;
Readonly my $SPACE       => q{ };
Readonly my $UNDERSCORE  => q{_};

### User Options
GetOptions(

    # Display help text
    q{h|help} => sub {
        done_testing();
        say {*STDOUT} $USAGE;
        exit;
    },

    # Display version number
    q{v|version} => sub {
        done_testing();
        say {*STDOUT} $NEWLINE
          . basename($PROGRAM_NAME)
          . $SPACE
          . $VERSION
          . $NEWLINE;
        exit;
    },
    q{vb|verbose} => $VERBOSE,
  )
  or (
    done_testing(),
    help(
        {
            USAGE     => $USAGE,
            exit_code => 1,
        }
    )
  );

BEGIN {

### Check all internal dependency modules and imports
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

## Modules
    my @modules = (q{MIP::QC::Record});

  MODULE:
    for my $module (@modules) {
        require_ok($module) or BAIL_OUT q{Cannot load} . $SPACE . $module;
    }
}

use MIP::QC::Record qw{ add_infile_info };

diag(   q{Test add_infile_info from Record.pm v}
      . $MIP::QC::Record::VERSION
      . $COMMA
      . $SPACE . q{Perl}
      . $SPACE
      . $PERL_VERSION
      . $SPACE
      . $EXECUTABLE_NAME );

## Create temp logger
my $test_dir = File::Temp->newdir();
my $test_log_path = catfile( $test_dir, q{test.log} );

## Creates log object
my $log = initiate_logger(
    {
        file_path => $test_log_path,
        log_name  => q{TEST},
    }
);

## Set file info parameters

# Interleaved files
my $is_interleaved;

my $read_length = $READ_LENGTH;
my $sample_id   = q{sample-1};

my %active_parameter = ( family_id => q{Adams}, );
my %file_info;
my %infile = ( $sample_id => [qw{ file-1 }], );
my %infile_both_strands_prefix;
my %infile_info = (
    date             => q{150703},
    direction        => 1,
    flowcell         => q{Undetermined-flow-rider},
    index            => q{ATCG},
    infile_sample_id => $sample_id,
    lane             => 1,
);
my %infile_lane_prefix;
my %lane;
my %sample_info;

## Define file formats
my ( $mip_file_format, $mip_file_format_with_direction,
    $original_file_name_prefix, $run_barcode )
  = _file_name_formats(
    {
        date      => $infile_info{date},
        direction => $infile_info{direction},
        flowcell  => $infile_info{flowcell},
        index     => $infile_info{index},
        lane      => $infile_info{lane},
        sample_id => $sample_id,
    }
  );

my $parsed_date = Time::Piece->strptime( $infile_info{date}, q{%y%m%d} );
$parsed_date = $parsed_date->ymd;

## Given single file, when Undetermined in flowcell name
SAMPLE_ID:
for my $sample_id ( keys %infile ) {

    my $lane_tracker = 0;

  INFILE:
    while ( my ( $file_index, $file_name ) = each @{ $infile{$sample_id} } ) {

        add_infile_info(
            {
                active_parameter_href           => \%active_parameter,
                date                            => $infile_info{date},
                direction                       => $infile_info{direction},
                file_index                      => $file_index,
                file_info_href                  => \%file_info,
                flowcell                        => $infile_info{flowcell},
                index                           => $infile_info{index},
                infile_both_strands_prefix_href => \%infile_both_strands_prefix,
                infile_href                     => \%infile,
                infile_lane_prefix_href         => \%infile_lane_prefix,
                is_interleaved                  => $is_interleaved,
                lane                            => $infile_info{lane},
                lane_href                       => \%lane,
                lane_tracker                    => $lane_tracker,
                read_length                     => $read_length,
                sample_id        => $infile_info{infile_sample_id},
                sample_info_href => \%sample_info,
            }
        );
    }
}

## Collect what to expect in one hash
my %expected_result = (
    lane => { $sample_id => [1], },
    infile_both_strands_prefix =>
      { $sample_id => [ $mip_file_format_with_direction, ], },
    infile_lane_prefix => {
        $sample_id => [ $mip_file_format, ],
    },
    sample_info => {
        sample => {
            $sample_id => {
                file => {
                    $mip_file_format => {
                        sequence_run_type   => q{single-end},
                        sequence_length     => $READ_LENGTH,
                        interleaved         => undef,
                        read_direction_file => {
                            $mip_file_format_with_direction => {
                                date               => $parsed_date,
                                original_file_name => q{file-1},
                                original_file_name_prefix =>
                                  $original_file_name_prefix,
                                read_direction => $infile_info{direction},
                                lane           => $infile_info{lane},
                                flowcell       => $infile_info{flowcell},
                                sample_barcode => $infile_info{index},
                                run_barcode    => $run_barcode,
                            },
                        },
                    },
                },
            },
        },
    },
);

## Then return true for detecting Undetermined in flowcell name
ok( $file_info{undetermined_in_file_name}{$mip_file_format},
    q{Tracked undetermined in file name} );

## Then add the lane info
is_deeply(
    \%lane,
    \%{ $expected_result{lane} },
    q{Added lane info for single-end read}
);

## Then add the infile lane prefix
is_deeply(
    \%infile_lane_prefix,
    \%{ $expected_result{infile_lane_prefix} },
    q{Added MIP file format for single-end read}
);

## Then add the infile both strands prefix (i.e. including read direction)
is_deeply(
    \%infile_both_strands_prefix,
    \%{ $expected_result{infile_both_strands_prefix} },
    q{Added MIP file format with direction for single-end read}
);

## Then add single-end read info from file name
is_deeply(
    \%sample_info,
    \%{ $expected_result{sample_info} },
    q{Added sample info for single-end read}
);

## Clear results
%expected_result = ();

## Given another file to mimic read direction 2
%infile = ( $sample_id => [qw{ file-1 file-2}], );

SAMPLE_ID:
for my $sample_id ( keys %infile ) {

    my $lane_tracker = 0;

  INFILE:
    while ( my ( $file_index, $file_name ) = each @{ $infile{$sample_id} } ) {

        if ( $file_index == 1 ) {

            ## Update infile info to mimic second read direction
            $infile_info{direction} = 2;

            ## Rdefine file format for new direction
            (
                $mip_file_format, $mip_file_format_with_direction,
                $original_file_name_prefix, $run_barcode
              )
              = _file_name_formats(
                {
                    date      => $infile_info{date},
                    direction => $infile_info{direction},
                    flowcell  => $infile_info{flowcell},
                    index     => $infile_info{index},
                    lane      => $infile_info{lane},
                    sample_id => $sample_id,
                }
              );

            ## Add sencond read to infile lane prefix
            push @{ $expected_result{infile_lane_prefix}{$sample_id} },
              $mip_file_format;
        }

        ## Add the infile both strands prefix
        push @{ $expected_result{infile_both_strands_prefix}{$sample_id} },
          $mip_file_format_with_direction;

        add_infile_info(
            {
                active_parameter_href           => \%active_parameter,
                date                            => $infile_info{date},
                direction                       => $infile_info{direction},
                file_index                      => $file_index,
                file_info_href                  => \%file_info,
                flowcell                        => $infile_info{flowcell},
                index                           => $infile_info{index},
                infile_both_strands_prefix_href => \%infile_both_strands_prefix,
                infile_href                     => \%infile,
                infile_lane_prefix_href         => \%infile_lane_prefix,
                is_interleaved                  => $is_interleaved,
                lane                            => $infile_info{lane},
                lane_href                       => \%lane,
                lane_tracker                    => $lane_tracker,
                read_length                     => $read_length,
                sample_id        => $infile_info{infile_sample_id},
                sample_info_href => \%sample_info,
            }
        );

        ## Define what to expect
        my %direction_one_metric = (
            interleaved     => $is_interleaved,
            sequence_length => $READ_LENGTH,
        );

        ## Alias
        my $file_level_href =
          \%{ $expected_result{sample_info}{sample}{$sample_id}{file}
              {$mip_file_format} };
      INFO:
        while ( my ( $file_key, $file_value ) = each %direction_one_metric ) {

            $file_level_href->{$file_key} = $file_value;
        }

        my %both_directions_metric = (
            date                      => $parsed_date,
            flowcell                  => $infile_info{flowcell},
            lane                      => $infile_info{lane},
            original_file_name        => $infile{$sample_id}[$file_index],
            original_file_name_prefix => $original_file_name_prefix,
            read_direction            => $infile_info{direction},
            run_barcode               => $run_barcode,
            sample_barcode            => $infile_info{index},
        );

        ## Alias
        my $direction_level_href =
          \%{ $expected_result{sample_info}{sample}{$sample_id}{file}
              {$mip_file_format}{read_direction_file}
              {$mip_file_format_with_direction} };

      INFO:
        while ( my ( $file_key, $file_value ) = each %both_directions_metric ) {

            $direction_level_href->{$file_key} = $file_value;
        }
    }
}

$expected_result{lane}{$sample_id} = [ 1, 1 ];

$expected_result{sample_info}{sample}{$sample_id}{file}{$mip_file_format}
  {sequence_run_type} = q{paired-end};

is_deeply(
    \%lane,
    \%{ $expected_result{lane} },
    q{Added lane info for paired-end read}
);

is_deeply(
    \%infile_lane_prefix,
    \%{ $expected_result{infile_lane_prefix} },
    q{Added MIP file format for paired-end read}
);
is_deeply(
    \%infile_both_strands_prefix,
    \%{ $expected_result{infile_both_strands_prefix} },
    q{Added MIP file format with direction for paired-end read}
);

is_deeply(
    \%sample_info,
    \%{ $expected_result{sample_info} },
    q{Added sample info for paired-end read}
);

done_testing();

######################
####SubRoutines#######
######################

sub build_usage {

## Function  : Build the USAGE instructions
## Returns   :
## Arguments : $program_name => Name of the script

    my ($arg_href) = @_;

    ## Default(s)
    my $program_name;

    my $tmpl = {
        program_name => {
            default     => basename($PROGRAM_NAME),
            store       => \$program_name,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    return <<"END_USAGE";
 $program_name [options]
    -vb/--verbose Verbose
    -h/--help     Display this help message
    -v/--version  Display version
END_USAGE
}

sub _file_name_formats {

## Function : Define format using information derived from infile name.
## Returns  : $mip_file_format, $mip_file_format_with_direction, $original_file_name_prefix, $run_barcode;
## Arguments: $date                            => Flow-cell sequencing date
##          : $direction                       => Sequencing read direction
##          : $flowcell                        => Flow-cell id
##          : $index                           => The DNA library preparation molecular barcode
##          : $lane                            => Flow-cell lane
##          : $sample_id                       => Sample id

    my ($arg_href) = @_;

    ## Flatten argument(s)
    my $date;
    my $direction;
    my $flowcell;
    my $index;
    my $lane;
    my $sample_id;

    my $tmpl = {
        date =>
          { defined => 1, required => 1, store => \$date, strict_type => 1, },
        direction => {
            allow       => [ 1, 2 ],
            defined     => 1,
            required    => 1,
            store       => \$direction,
            strict_type => 1,
        },
        flowcell => {
            defined     => 1,
            required    => 1,
            store       => \$flowcell,
            strict_type => 1,
        },
        index =>
          { defined => 1, required => 1, store => \$index, strict_type => 1, },
        lane => {
            allow       => qr/ ^\d+$ /xsm,
            defined     => 1,
            required    => 1,
            store       => \$lane,
            strict_type => 1,
        },
        sample_id => {
            defined     => 1,
            required    => 1,
            store       => \$sample_id,
            strict_type => 1,
        },
    };

    check( $tmpl, $arg_href, 1 ) or croak q{Could not parse arguments!};

    my $mip_file_format =
        $sample_id
      . $DOT
      . $date
      . $UNDERSCORE
      . $flowcell
      . $UNDERSCORE
      . $index
      . $DOT . q{lane}
      . $lane;

    my $mip_file_format_with_direction =
      $mip_file_format . $UNDERSCORE . $direction;

    my $original_file_name_prefix =
        $lane
      . $UNDERSCORE
      . $date
      . $UNDERSCORE
      . $flowcell
      . $UNDERSCORE
      . $sample_id
      . $UNDERSCORE
      . $index
      . $UNDERSCORE
      . $direction;

    my $run_barcode =
        $date
      . $UNDERSCORE
      . $flowcell
      . $UNDERSCORE
      . $lane
      . $UNDERSCORE
      . $index;
    return $mip_file_format, $mip_file_format_with_direction,
      $original_file_name_prefix, $run_barcode;
}
