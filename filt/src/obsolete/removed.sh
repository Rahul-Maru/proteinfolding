# OBSOLETE
# removes binding sites that are not present in the PDB after an rysnc
#   but were not removed from the binding-sites dir

shopt -s nullglob

for d in binding-sites/*; do
	echo $d

	for f in "$d"/*; do
		nam=$( basename $f )
		id=${nam:3:4}
		if ! grep -iq $id pdb.txt; then
			rm "$f"
			echo "removed $f"
		fi
	 done
done
