for dir in pdb/*; do
	d=$(basename "$dir")
	echo "$d"
	mkdir "binding-sites/$d"
done
for f in binding-sites/*; do
	if [ -f "$f" ]; then
		name=$(basename "$f")
		mid=${name:4:2}
		mv "$f" "binding-sites/$mid"
		echo "$name moved to $mid/"
	fi
done

shutdown

