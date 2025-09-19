inf="xa$1"
echo "reading from $inf"

mkdr=:

if [ "$mkdr" = true ]; then
	for dir in pdb/*; do
		d=$(basename "$dir")
		echo "$d"
		mkdir "binding-sites/$d"
	done
fi


while read f; do
	f="binding-sites/$f"
	if [ -f "$f" ]; then
		name=$(basename "$f")
		mid=${name:4:2}
		mv "$f" "binding-sites/$mid"
	fi
done < "$inf"

c=$(ls -p "binding-sites" | grep -v "/" | wc -w)
echo "$c"
if [ $c -eq 0 ]; then
	shutdown
fi


