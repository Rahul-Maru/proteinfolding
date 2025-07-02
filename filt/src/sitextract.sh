nam=$(basename $1)
mid=${nam:1:2}
echo "$mid - $nam"
if [ -f "../dat/struct/pdb/$mid/pdb$1.ent" ]; then
    echo "pdb$1.ent found"
fi
./../dat/struct/sitextract -o test "../dat/struct/pdb/$mid/pdb$1.ent"
gedit "../dat/struct/pdb/$mid/pdb$1.ent"