cd tools/MAPP-3D/MultipleSiteAlignment/

for i in {1..19}; do
	for j in {1..10}; do
		echo -n "$(calc $i/20); $j; "
		python Analyse.py align_output.txt $(calc $i/20) $j
	done
done
echo

