#!/bin/bash

shopt -s nullglob

find ./filtered-binding-sites/ -name "*.pdb" -type f -delete

for d in binding-sites/*; do
	echo "d - $d"
	for f in $d/*; do
		nam=$( basename $f )	
		ligArr=(${nam//_/ })
		lig=${ligArr[3]}
		lig=${lig%.*}

		len=$( awk '{print substr($0, 22,5)}' $f | sort -n | uniq | wc -l )
		if (( len >= 4 )) && ! grep -qw $lig unwanted-ligs; then
			mid=${nam:4:2}
			cp $f filtered-binding-sites/$mid/$nam
		fi
	done 
done

