#!/usr/bin/env/bash

usage="$(basename "$0") [-e <environment_name> -p <conda_path>]

where:
    -h  Show this help text
    -e  Conda environment name
    -p  Conda path"

unset OPTARG
unset OPTIND

ENV_NAME="mip7_rd-dna"
CONDA_PATH="$HOME/miniconda3"

while getopts ':he:p:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    e) ENV_NAME=$OPTARG
       ;;
    p) CONDA_PATH=$OPTARG
       ;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
done
shift $((OPTIND-1))

conda create --name "$ENV_NAME" --yes  -c conda-forge libgcc-ng gxx_linux-64

conda install --name "$ENV_NAME" --yes -c conda-forge perl=5.26 perl-app-cpanminus

conda install --name "$ENV_NAME" --yes -c bioconda perl-log-log4perl perl-moosex-app

## Source conda
source "$CONDA_PATH"/etc/profile.d/conda.sh

conda activate "$ENV_NAME"

cpanm --installdeps .

conda deactivate

