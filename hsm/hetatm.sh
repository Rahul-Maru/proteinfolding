#!/bin/bash

indir=$1
outdir=$2

for f in $indir/*; do
	ext=${f##*.}
	if [[ $ext == "pdb" ]]; then
		echo $f
		nam=$(basename $f)
		awk '
		function get_lig(l) {
			return substr(l, 18,3)
		}

		!/^ATOM/ { print $0 }

		/^ATOM/ {
			lin = $0
			lig = get_lig($0)
			if (lig == "HSM") {
				sub(/^ATOM  /, "HETATM")
				print
			} else {
				print $0
			}
		}
		' $f > "$outdir/$nam"
	fi
done
