inf=$1
outf=$2

grep -E "^ATOM|^HETATM\s+.*(ZN|FE|MG)|^TER|^END" $inf > $outf site.pdb