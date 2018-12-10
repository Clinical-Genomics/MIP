#!/usr/bin/env/bash
set -e
set -u

mkdir --parents /Users/henrik.stranneheim/git/MIP/.test_create_housekeeping_function 

finish() {

	local directory="$1"
	## Perform exit housekeeping
	rm --recursive --force "$directory" 

}

## Enable trap for signal(s) EXIT TERM INT
trap '$(finish /Users/henrik.stranneheim/git/MIP/.test_create_housekeeping_function)' EXIT TERM INT 

rm /Users/henrik.stranneheim/git/MIP/test_create_housekeeping_function.sh 