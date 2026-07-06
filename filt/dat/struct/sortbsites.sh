#!/bin/bash

# sorts binding-sites into subdirectories as in pdb/
#   defined by the middle 2 chars of the pdb id
#   takes in as input one of the outputs of 'split'
#	pass 'sd' as second input to shut down after all bsites are moved

# next step: filterbsites.sh

inf="$1"
indir="binding-sites"

if ! [ -f $inf ]; then
	echo "input the slice output file to be read (or bsites.txt)"
	echo "eg. ./sortbsites.sh xaa"
	exit 1
fi

echo "reading from $inf"

# if the dir has no subdirs then make them
subdircount=$(find "$indir" -maxdepth 1 -type d | wc -l)
if [[ "$subdircount" -eq 1 ]] ; then
	echo making dirs
	for dir in pdb/*; do
		d=$(basename "$dir")
		echo "$d"
		mkdir "$indir/$d"
	done
fi


while read f; do
	f="$indir/$f"
	if [ -f "$f" ]; then
		name=$(basename "$f")
		mid=${name:4:2}

		mv "$f" "$indir/$mid"
	fi
done < "$inf"

c=$(ls -p "$indir" | grep -v "/" | wc -w)
echo "$inf sorted"
echo "$c sites left"

if [ $c -eq 0 ]; then
	[[ $2 == "sd" ]] && shutdown
fi


