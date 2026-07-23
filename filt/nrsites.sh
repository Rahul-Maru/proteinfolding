#!/bin/bash

# MASTER SCRIPT TO CURATE DATABASE OF NON-REDUNDANT BINDING SITE STRUCTURES FROM THE PDB
# WRITTEN BY RAHUL MARU, AS PART OF THE NSc LAB IN THE DEPT. OF BIOCHEMISTRY, IISc
#
# SCRIPTS CAN ALSO BE RUN INDIVIDUALLY
#
# USAGE:
#	script must be run from the bio/filt/ directory
#   . nrsites.sh            run the whole pipeline
#   . nrsites.sh <step>     start from <step> and run everything after it

shopt -s nullglob

root=$(pwd)
n=4   # number of parallel workers; raise on a bigger machine

timed() {
	t0=$SECONDS
	"$@"
	t=$(($SECONDS - $t0))

	hrs=$(( t / 3600))
	mns=$(( t / 60 % 60 ))
	scs=$(( t % 60 ))

	name=$1
	[[ $1 == "python" ]] && name=$2
	echo "$name done in ${hrs}h, ${mns}m, ${scs}s"
	echo
}

move_old() {
	cd "$root/dat/struct" || return 1
	local archive="old3/$(date +%Y%m%d-%H%M%S)"
	mkdir -p "$archive"
	echo "archiving existing directories to filt/dat/struct/$archive/"
	for d in */; do
		[[ "$d" == old3/ ]] && continue
		mv "$d" "$archive/"
	done
	mv NRsites.tar.gz "$archive/"
	printf "done\n\n"
}

cluster_file_download() {
	cd "$root/dat" || return 1
	echo "downloading cluster-file"
	timed curl -sO https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-70.txt
	printf "done\n\n"
	echo "cluster-file downloaded on $(date +%d/%m/%Y) (dd/mm/yy)" > ../README.md
}

mkdirs() {
	cd "$root/dat/struct" || return 1
	mkdir -p binding-sites filtered-binding-sites representative-binding-sites
	printf "made directories\n\n"
}

pdb() {
	cd "$root/dat/struct" || return 1
	echo "downloading pdb via rsync"
	timed "$root/src/rsyncPDB.sh"

	echo "PDB files downloaded via rsync on $(date +%d/%m/%Y)" >> ../../README.md
}

unzip() {
	cd "$root/dat/struct" || return 1
	echo "unzipping files"
	timed "$root/src/extractzip.sh"
}

extract() {
	echo "extracting molecular entity numbers into $root/dat/struct/ents2.txt"
	cd "$root" || return 1
	timed python src/extractents.py

	cd "$root/dat/struct" || return 1
	echo "extracting binding sites into $(pwd)/binding-sites"
	timed "$root/src/bsites.sh"

	echo "listing binding sites into $(pwd)/bsites.txt"
	ls -R binding-sites | grep pdb > bsites.txt
}

sortbsites() {
	cd "$root/dat/struct" || return 1

	if [[ ! -f bsites.txt ]]; then
		echo "bsites.txt not found — run the 'extract' step first"
		return 1
	fi

	echo "pre-creating subdirectories in $(pwd)/binding-sites/"
	"$root/src/mksubdir.sh" binding-sites

	echo "sorting binding sites with $n parallel workers"
	split -n "l/$n" bsites.txt sortslice.

	# wrapper for timing purposes
	sorting_bsites() {
		for slice in sortslice.*; do
			"$root/src/sortbsites.sh" "$slice" &
		done
		wait
	}

	timed sorting_bsites
	rm -f sortslice.*
}

filter() {
	cd "$root/dat/struct" || return 1
	echo "creating subdirectories for filtering"
	"$root/src/mksubdir.sh" filtered-binding-sites
	printf "done\n\n"

	echo "filtering binding sites"
	cd "$root" || return 1
	timed python src/filterbsites.py

	cd "$root/dat/struct" || return 1
	ls -R filtered-binding-sites | grep pdb > f_bsites.txt
}

clusterize() {
	cd "$root" || return 1
	timed python src/clusterize.py
}

reprs() {
	cd "$root" || return 1
	timed python src/choose_reprs.py
}

move_reprs() {
	echo "converting output .json file ($root/dat/reprs.json)
		into .txt ($root/dat/struct/cluster_reprs.txt)"
	cd "$root" || return 1
	timed python src/reprsjson_to_txt.py

	cd "$root/dat/struct" || return 1
	echo "moving representative binding sites to $(pwd)/representative-binding-sites"
	timed "$root/src/move_reprs.sh"
}

move_lig() {
	echo "copying ligands into $root/dat/struct/ligand-inclsv-binding-sites"
	cd "$root" || return 1
	timed python src/lig.py
}

compress() {
	cd "$root/dat/struct" || return 1
	echo "compressing ligand-inclsv-binding-sites/"
	timed tar czf NRsites.tar.gz ligand-inclsv-binding-sites
	mkdir -p "$root/out"
	mv NRsites.tar.gz ../../out/
	printf "\nNRsites COMPLETE. THE FINAL ARCHIVE CAN BE FOUND AT filt/out/NRsites.tar.gz. THE UNARCHIVED FOLDER CAN BE FOUND AT filt/dat/struct/ligand-inclsv-binding-sites. TO VERIFY OUTPUT, RUN THE SCRIPTS IN filt/src/analysis, BEARING IN MIND THAT SOME ARE OBSOLETE.\n"
}

# input validation

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
	echo "script must be run via: . nrsites.sh  NOT:  ./nrsites.sh"
	exit 1
fi

if [[ ! -d "$root/dat/struct" ]]; then
	echo "dat/struct directory not found"
	mkdir -p "$root/dat/struct"
fi


# driver loop -- calls all the functions

steps=(move_old cluster_file_download mkdirs pdb unzip extract sortbsites filter clusterize reprs move_reprs move_lig compress)

start="${1:-${steps[0]}}"

if [[ ! " ${steps[*]} " == *" $start "* ]]; then
	echo "unknown start step: $start"
	echo "valid steps: ${steps[*]}"
	return 1
fi

seen=0
for step in "${steps[@]}"; do
	[[ "$step" == "$start" ]] && seen=1
	(( seen )) && "$step"
done

cd "$root" || return 1

