INDR=""

for f in $INDR/*; do
	nam=$( basename $f )	
	ligArr=(${nam//_/ })
	n=${ligArr[2]}
	og_nam="${ligArr[0]}.pdb"
	echo $n = $og_nam

	# awk '
	# resn = substr($0, 23, 4)
	# if 
	# '

done