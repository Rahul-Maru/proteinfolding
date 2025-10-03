for d in binding-sites/*; do
	for f in "$d"/*; do
		if [ -f "$f" ] && [ ! -s "$f" ]; then
			echo "removed $f"
			rm $f
		fi
	done
done
