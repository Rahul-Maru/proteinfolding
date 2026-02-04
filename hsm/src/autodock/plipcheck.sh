cd ../../outs/autodock/poses

for f in *; do
    mkdir ../plip/$f
    ./../../../tools/plip.simg -f $f -o ../plip/$f/ -vy
done