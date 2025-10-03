odir="binding-sites"
cutoff=4.5
for d in pdb/*; do
	echo "$d"
	for f in "$d"/*; do
		./sitextract -o "$odir" -d "$cutoff" "$f"
	done
done

