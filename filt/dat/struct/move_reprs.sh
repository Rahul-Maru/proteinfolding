# copies chosen representatives into new dir

# empties directory first
find "representative-binding-sites" -type f -delete

while IFS= read -r repr; do

	mid=${repr:4:2}

	cp "filtered-binding-sites/$mid/$repr" "representative-binding-sites/$repr"

done < cluster_reprs.txt
