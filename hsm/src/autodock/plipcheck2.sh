cd ../../outs/autodock/poses/
echo "running from $(pwd)"

for lig in *.pdbqt; do
    nam=${lig::-6}
    echo "----processing $nam----"

    recnam="${nam:4:2}/${nam:0:7}.ent"
    rec="../../../../filt/dat/struct/pdb/$recnam"

    d="../plip/$nam"
    tmplig="$d/tmp_lig_$nam.pdb"
    tmplig2="$d/tmp2_lig_$nam.pdb"

    tmprec="$d/tmp_rec_$nam.pdb"
    tmpf="$d/tmp_$nam.pdb"
    tmpf_fix="$d/temp_${nam}_fixed.pdb"

    if [ -d $d ]; then
        rm $d/*
        echo "cleared $d"
    else
        mkdir "$d/"
        echo "made $d"
    fi

    # fix receptor
    python -c "from bio import fixpdb; fixpdb('$rec', '$tmprec')"
    grep -E "ATOM" $tmprec > $tmpf
    echo "fixpdb done"

    obabel -i pdbqt $lig -o pdb -O $tmplig
    echo "obabel done"

    grep "HETATM" $tmplig >> $tmpf

    pdb_fixinsert $tmpf | pdb_reatom | pdb_reres > $tmpf_fix

    echo "fixed insertions and renumbered atoms"

    ./../../../tools/plip2.simg -f $tmpf_fix -o $d/ -vyt
    echo "plip done"

    rm $tmplig $tmprec $tmpf

    echo ----done with $nam----
    echo
done
