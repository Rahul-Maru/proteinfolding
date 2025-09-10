#!/bin/bash

INDR="representative-binding-sites"
OUTDIR="ligand-inclsv-binding-sites"

if ! [ -d $OUTDIR ]; then
	echo $OUTDIR
	mkdir $OUTDIR
fi

for f in $INDR/*; do
	nam=$( basename $f )
	ligArr=(${nam//_/ })
	n=${ligArr[2]}
	og_nam="${ligArr[0]}.ent"
	subdir=${og_nam:4:2}
	# echo $nam: $n = ./pdb/$subdir/$og_nam

	found=0
	while IFS= read -r line; do
		rec=${line:0:6}
		if [[ $rec == "HETATM" ]]; then

			num=${line:22:4}
			num="${num//[[:space:]]/}"
			num=$(printf "%04d" "$num")
			if [[ "$num" == "$n" ]]; then
				if [[ $found == 0 ]]; then
					found=1
					if ! [ -f $OUTDIR/$nam ]; then
						cp $f $OUTDIR/$nam
					fi
				fi
				echo $line >> $OUTDIR/$nam
			fi
		fi
	done < "pdb/$subdir/$og_nam"
	if [[ $found == 0 ]]; then
		echo "not found $nam"
	fi
done
