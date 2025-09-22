#!/bin/bash

inf="xa$1"
echo "reading from $inf"

mkdr=""

if [[ -n "$mkdr" ]]; then
	echo hi
	for dir in pdb/*; do
		d=$(basename "$dir")
		echo "$d"
		mkdir "binding-sites/$d"
	done
fi


while read f; do
	f="binding-sites/$f"
	if [ -f "$f" ]; then
		name=$(basename "$f")
		mid=${name:4:2}
		if ! [ -d "binding-sites/$mid" ]; then
			echo oh noes $mid
			exit
		fi
		mv "$f" "binding-sites/$mid"
	fi
done < "$inf"

c=$(ls -p "binding-sites" | grep -v "/" | wc -w)
echo "$c"
if [ $c -eq 0 ]; then
	shutdown
fi


