while IFS= read -r repr; do
	mid=${repr:4:2}

	cp "binding-sites/$mid/$repr" "representative-binding-sites/$repr"

done < cluster_reprs.txt
