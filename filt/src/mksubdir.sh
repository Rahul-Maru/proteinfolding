# creates subdirectories in the base directory given as the first argument

echo "making dirs in $1/"

for dir in pdb/*; do
	d=$(basename "$dir")
	echo "$d"
	mkdir -p "$1/$d"
done
