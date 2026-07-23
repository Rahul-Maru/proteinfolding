#!/bin/bash

# removes empty binding sites

for d in binding-sites/*; do
	for f in "$d"/*; do
		if [ -f "$f" ] && [ ! -s "$f" ]; then
			echo "removed" $(basename $f)
			rm $f
		fi
	done
done
