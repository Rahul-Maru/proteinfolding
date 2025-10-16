indr="tools/MAPP-3D/MultipleSiteAlignment/flappsites"

if ! [ -d $indr ]; then
    mkdir $indr
    echo $indr
else
    find $indr -type f -delete
fi

while read f; do
    echo $f
    cp "tools/FLAPP/sites/$f" "$indr/"
done < outs/FLAPP/alignlist.txt

for s in tools/FLAPP/sites/*HSM*; do
    echo $s
    cp "$s" "$indr/"
done