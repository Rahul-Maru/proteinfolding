odir="binding-sites"
cutoff=4.5
for d in pdb/*; do
	echo "$d"
	for f in "$d"/*; do
		nam=$(basename $f)
		nam=${nam:0:7}
		mid=${nam:4:2}
		if [[ $(grep "$nam" "bsites.txt" | wc -w) == 0 ]]; then
			echo "uh $nam"
			./sitextract -o "$odir/$mid" -d "$cutoff" "$f"
		fi
	done
done
