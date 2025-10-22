# creates subdirectories

echo "making dirs"

for dir in pdb/*; do
	d=$(basename "$dir")
	echo "$d"
	mkdir "$1/$d"
done
