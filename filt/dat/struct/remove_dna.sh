#!/bin/bash

shopt -s nullglob

filter_file() {
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

	!/^ATOM/ { print $0 }

	/^ATOM/ {
		res = substr($0, 18, 3)
		if (res in aacodes) {
			print $0
		}
	}' "$f"
}

bsites_base='/home/aibio/Desktop/bio/filt/dat/struct/binding-sites'
out_base='/home/aibio/Desktop/bio/filt/dat/struct/filtered-binding-sites'

# find ./$outbase/ -name "*.pdb" -type f -delete


if [ ! -d "$out_base" ]; then
	mkdir -p "$out_base"
fi

for d in "${bsites_base}/"*; do
	if [ -d "$d" ]; then
		echo "d - $d"
		dbase=$(basename "$d")
		outd="${out_base}/${dbase}"
		if [ ! -d "$outd" ]; then
			mkdir -p "$outd"
		fi
		for f in "${d}/"*.pdb; do
			fbase=$(basename "$f")
			ligArr=(${fbase//_/ })
			lig=${ligArr[3]}
			lig=${lig%.*}
			len=$(filter_file "$f" | awk '/^ATOM/ {print substr($0, 22,5)}' | sort -n | uniq | wc -l )
			outf="${outd}/${fbase}"
			if (( len >= 4 )) && ! grep -qw $lig unwanted-ligs.txt; then
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
