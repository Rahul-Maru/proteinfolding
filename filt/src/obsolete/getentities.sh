#!/bin/bash

# OBSOLETE (superseded by extractents.py)

shopt -s nullglob

outf="ents.txt"
errf="noents.txt"

# clears output and error files
> $outf
> $errf

while IFS= read -r pdb; do
	i=0
	nam=${pdb:3:4}
	mid=${nam:1:2}
	#echo ----$nam----
	while IFS= read -r line; do
		if mat=$(echo $line | grep -oP "^COMPND\s+\d*\s*MOL_ID:\s*(\d+)"); then
			#echo found: $mat

			# gets the entity number and writes it to outf
			ent=$(echo $mat | cut -d : -f 2 | grep -oP "\d*")
			#echo $i-${nam}_${ent}
			echo ${nam}_${ent} >> $outf
			i=$((i+1))
			#echo found in: $line
		elif echo $line | grep -q "^SOURCE"; then
			#echo src: $line
			if [[ $i == 0 ]]; then
				echo $nam has no entities
				echo $nam >> $errf
			fi
			break
		fi
	done < pdb/$mid/$pdb
done < pdb.txt

