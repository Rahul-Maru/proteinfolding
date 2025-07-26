#!/bin/bash

for d in binding-sites/*; do
	echo "d - $d"
	for f in $d/*; do
		nam=$( basename $f )	
		echo "n - $nam"

		ligArr=(${f//_/ })
		lig=${ligArr[3]}
		lig=${lig%.*}

		if grep -qw $lig unwanted-ligs; then
			echo "l - $lig"
		fi
		echo "done"
	done 
done
