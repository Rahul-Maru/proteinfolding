#!/bin/bash

# OBSOLETE (superseded by filterbsites.sh)

# binding sites are copied into filtered-binding-sites if:
#	1. they do not bind to an unwanted ligand
#   AND
#	2. they have at least 4 residues

shopt -s nullglob

find ./filtered-binding-sites/ -name "*.pdb" -type f -delete

for d in binding-sites/*; do
	echo "d - $d"
	for f in $d/*; do
		nam=$( basename $f )	

		# ligand
		ligArr=(${nam//_/ })
		lig=${ligArr[3]}
		lig=${lig%.*}

		# number of residues
		len=$( awk '{print substr($0, 22,5)}' $f | sort -n | uniq | wc -l )

		if (( len >= 4 )) && ! grep -qw $lig unwanted-ligs.txt; then
			mid=${nam:4:2}
			cp $f filtered-binding-sites/$mid/$nam
		fi
	done 
done

