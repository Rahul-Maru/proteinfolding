#!/bin/bash

# OBSOLETE (superseded by filterbsites.py)

# filters out unwanted binding sites.
# binding sites are removed if:
#	1. they bind to an unwanted ligand
#   OR
#	2. they have fewer than 4 amino acid residues

# TODO speed this up

shopt -s nullglob

if [[ $(basename $(pwd)) != "struct" ]]; then
	echo "Error: must be run from the struct/ directory"
	exit 1
fi

if [ ! -d "$out_base" ]; then
	mkdir -p "$out_base"
fi

# if the dir has no subdirs then make them
subdircount=$(find "filtered-binding-sites" -maxdepth 1 -type d | wc -l)
if [[ "$subdircount" -eq 1 ]]; then
	for dir in pdb/*; do
		d=$(basename "$dir")
		echo "$d"
		mkdir "filtered-binding-sites/$d"
	done
fi

filter_file() {
	# function to filter non-AA residues
	local f="$1"
	awk '
	BEGIN {
		aacodes["ALA"] = 1
		aacodes["ARG"] = 1
		aacodes["ASN"] = 1
		aacodes["ASP"] = 1
		aacodes["CYS"] = 1
		aacodes["GLU"] = 1
		aacodes["GLN"] = 1
		aacodes["GLY"] = 1
		aacodes["HIS"] = 1
		aacodes["ILE"] = 1
		aacodes["LEU"] = 1
		aacodes["LYS"] = 1
		aacodes["MET"] = 1
		aacodes["PHE"] = 1
		aacodes["PRO"] = 1
		aacodes["SER"] = 1
		aacodes["THR"] = 1
		aacodes["TRP"] = 1
		aacodes["TYR"] = 1
		aacodes["VAL"] = 1
	}

	# there should be no non-ATOM records in input files
	!/^ATOM/ { print $0 }

	/^ATOM/ {
		# gets the res code
		res = substr($0, 18, 3)
		if (res in aacodes) {
			print $0
		}
	}' "$f"
}

bsites_base='binding-sites'
out_base='filtered-binding-sites'
unwanted_ligs='unwanted-ligs'

# find ./$outbase/ -name "*.pdb" -type f -delete

for d in "${bsites_base}/"*; do
	if [ -d "$d" ]; then
		echo "d - $d"
		dbase=$(basename "$d")
		outd="${out_base}/${dbase}"

		# make subdirectory
		if [ ! -d "$outd" ]; then
			mkdir -p "$outd"
		fi

		for f in "${d}/"*.pdb; do
			fbase=$(basename "$f")

			# get ligand
			ligArr=(${fbase//_/ })
			lig=${ligArr[3]}
			lig=${lig%.*}

			# get number of residues after removing non-AAs
			len=$(filter_file "$f" | awk '/^ATOM/ {print substr($0, 22,5)}' | sort -n | uniq | wc -l )

			outf="${outd}/${fbase}"
			if (( len >= 4 )) && ! grep -qw "$lig" "unwanted-ligs.txt"; then
				filter_file "$f" > "$outf"
			else
				if [ -f "$outf" ]; then
					echo "removing $outf"
					rm "$outf"
				fi
            fi
		done
	fi
done
