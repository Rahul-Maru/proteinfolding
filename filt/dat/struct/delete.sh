for d in binding-sites/*; do
	for f in "$d/*"; do
		echo $f
	done
done
