#!/bin/bash

# extract all binding sites from pdb/ into odir

# next step:
# * ls -R binding-sites | grep pdb > bsites.txt
# * split -n l/[4] bsites.txt (4 can be changed)
# * ./sortbsites.sh [inf] (where inf is each output of split)

shopt -s nullglob

# sitextract lives alongside this script in src/; find it by script path, not cwd
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

odir="binding-sites"
cutoff=4.5
for d in pdb/*; do
	echo "$d"
	for f in "$d"/*; do
		"$script_dir/sitextract" -o "$odir" -d "$cutoff" "$f"
	done
done

