#!/bin/bash

for f in cif/*; do
	nam=$( basename $f )
	id=${nam:0:4}
	if ! [ -f "pdbs2/$id.pdb" ]; then
		echo $id
	else
		rm $f	
	fi
done
