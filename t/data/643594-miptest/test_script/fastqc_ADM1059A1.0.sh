#!/bin/bash --login
set -e
set -o pipefail

#SBATCH --account=cust00X
#SBATCH --ntasks=1
#SBATCH --time=00:00:00
#SBATCH --qos=low
#SBATCH --job-name=fastqc_ADM1059A1
#SBATCH --error=643594-miptest/analysis/ADM1059A1/fastqc/info/fastqc_ADM1059A1.0.stderr.txt
#SBATCH --output=643594-miptest/analysis/ADM1059A1/fastqc/info/fastqc_ADM1059A1.0.stdout.txt

readonly PROGNAME=$(basename "$0")

echo "Running on: $(hostname)" 

## Create temporary directory
readonly TEMP_DIRECTORY=/scratch/"$SLURM_JOB_ID"
mkdir --parents "$TEMP_DIRECTORY" 

finish() {

	local directory="$1"
	## Perform exit housekeeping
	rm --recursive --force "$directory" 

}

## Enable trap for signal(s) EXIT TERM INT
trap '$(finish "$TEMP_DIRECTORY")' EXIT TERM INT 


## Enable trap for signal(s) DEBUG
trap 'previous_command="$BASH_COMMAND"' DEBUG 

error() {

	local program="$1"

	local return_code="$2"

	## Display error message and exit
	echo "${program}: ${return_code}: Unknown Error - ExitCode=$return_code" 1>&2
	exit 1
}

## Enable trap for signal(s) ERR
trap '$(error "$previous_command" "$?")' ERR 


