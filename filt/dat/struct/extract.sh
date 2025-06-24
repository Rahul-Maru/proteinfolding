for d in pdb/*; do
	if [ -d "$d" ]; then
		echo "$d"
		for f in "$d"/*; do
			if [ $(echo "$f" | grep -c "ent[.]gz$") -eq 1 ]; then
				echo "hi"
				gunzip "$f"
			fi
		done
	fi
done
