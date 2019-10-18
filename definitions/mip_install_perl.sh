#!/usr/bin/env/bash
set -e

usage="$(basename "$0") [-e <environment_name> -p <conda_path> -x <existing_conda_env>]

where:
    -h  Show this help text
    -e  Conda environment name [string]
    -p  Conda path [string]
    -x  Install in already existing conda environment [flag]"

unset OPTARG
unset OPTIND

CONDA_PATH="$HOME/miniconda3"
ENV_NAME='mip7_rd-dna'
EXISTING_ENV='false'

while getopts ':he:p:x' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    e) ENV_NAME=$OPTARG
       ;;
    p) CONDA_PATH=$OPTARG
       ;;
    x) EXISTING_ENV='true'
       ;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
done
shift $((OPTIND-1))

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

## Create or install conda env
if [ "$EXISTING_ENV" = true ]
then
  conda install --name "$ENV_NAME" --yes  -c conda-forge libgcc-ng gxx_linux-64
else
  conda create --name "$ENV_NAME" --yes  -c conda-forge libgcc-ng gxx_linux-64
fi

conda install --name "$ENV_NAME" --yes -c conda-forge perl=5.26 perl-app-cpanminus

conda install --name "$ENV_NAME" --yes -c bioconda perl-log-log4perl perl-moosex-app

## Source conda
source "$CONDA_PATH"/etc/profile.d/conda.sh

conda activate "$ENV_NAME"

cpanm --installdeps .

conda deactivate

