#!/bin/bash

make_dirs() {
#	if [ -z "$2" ]; then
#		indir=""
#		loop="old"
#	else 
#		indir=$2
#		loop="old$indir"
#	fi

	target=${2:4}
	if ! [ -d $target ]; then
		echo making $target
		mkdir $target
	fi	

	for d in $2/*; do
		if [ -d $d ]; then
		#	nam=$(basename $d)
		#	echo "in loop: $indir - $d | $indir/$nam"

			j=$(( $1 + 1 ))
			for ((i=0;i<j;i++)); do echo -n "--"; done

			make_dirs $j "$d"
		#	echo onwards
		fi
	done
}

echo .
make_dirs 0 old
