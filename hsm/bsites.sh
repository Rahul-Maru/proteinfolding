odir="bsites-all"
cutoff=5.5
for f in pdbs/*; do
	echo "$f"
	nam=$(basename $f)
	nam=${nam:0:4}
	mid=${nam:1:2}
	mkdir "$odir/$mid"
	./sitextract -o "$odir/$mid" -d "$cutoff" "$f"
done
