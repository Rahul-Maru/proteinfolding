for d in binding-sites/*; do
	for f in "$d"/*; do
		# Test the negation
		if [ -f "$f" ] && [ ! -s "$f" ]; then
			echo "removed $f"
			rm $f
		fi
	done
done
