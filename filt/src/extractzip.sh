# unzips all pdb files, n at a time
# next step: bsites.sh

n=${1:-4}   # number of parallel workers (first arg, default 4)

find pdb -type f -name '*.ent.gz' -print0 | xargs -0 -P "$n" -n 1 gunzip
