#!/bin/bash

for f in cif/*; do
	nam=$(basename $f)
	nam=${nam:0:4}
	obabel -i cif $f -o pdb cif/ -O "cif/$nam.pdb"
done
