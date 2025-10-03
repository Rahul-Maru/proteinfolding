indir=$1

for f in $indir/*; do
	gunzip -dv $f
done

