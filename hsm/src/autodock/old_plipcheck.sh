# OBSOLETE. USE PLIPCHECK2

echo "OBSOLETE WARNING"

cd ../../outs/autodock/poses/

for f in *.pdbqt; do
    nam="${f::-6}"
    d="../plip/$nam"
    tmpf="$d/tmp_$nam.pdb"

    mkdir "$d/"
    obabel -i pdbqt $f -o pdb -O $tmpf
    ./../../../tools/plip.simg -f $tmpf -o $d/ -vyt
    rm $tmpf
done
