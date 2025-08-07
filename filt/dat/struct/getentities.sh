shopt -s nullglob

outf="ents.txt"

for d in pdb/*; do
	for f in "$d"/*; do
		nam=$( basename $f )
		pdb=${nam:3:4}
		if ! grep -iq $pdb "../clusters-by-entity-70.txt"; then
			echo $pdb - $nam
		fi
		#while IFS= read -r line; do
		#	break			
		#done < pdb.txt
	done
done

