# extract all binding sites from pdb/ into odir

# next step: ls -R binding-sites | grep pdb > bsites.txt
# * split -n 4 bsites.txt
# * sortbsites.sh

odir="binding-sites"
cutoff=4.5
for d in pdb/*; do
	echo "$d"
	for f in "$d"/*; do
		./sitextract -o "$odir" -d "$cutoff" "$f"
	done
done

