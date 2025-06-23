dir="bchains_seq"
dir2="bsites_final"
echo "hi"
i=$(ls "$dir" | wc -w)
j=$(ls "$dir2" | wc -w)
echo "dir1: $dir-$i"
echo "dir2: $dir2-$j"

for f in "$dir"/*; do
	fnam=$(basename "$f" | sed 's/\.pdb//')
	
	flg=""
	if ls "$dir2" | grep "$fnam"; then
		flg="okay"
	else
		flg="not okay"
	fi
	echo "e - $fnam - $flg"
done

