cd tools/MAPP-3D/MultipleSiteAlignment/

for i in {1..16}; do
	for j in {1..8}; do
		echo -n "$(calc $i/20), $j, "
		python Analyse.py align_output.txt $(calc $i/20) $j
	done
done
echo

