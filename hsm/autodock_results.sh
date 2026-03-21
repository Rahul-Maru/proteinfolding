for d in outs/autodock/results/*/; do
	if [ -f "$d/1.dlg" ]; then
		e=$(cat "$d/1.dlg" | grep RANKING | head -n 1 | awk '{print $4}')
		s=$(basename $d)
		echo "$s, $e"
	fi
done
