# COPIES CHOSEN REPRESENTATIVES INTO NEW DIR

# makes sure that relevant files and directories exist before emptying representative-binding-sites
if [[ ! -s cluster_reprs.txt ]]; then
	echo "Error: cluster_reprs.txt is missing or empty"; exit 1
fi

mkdir -p "representative-binding-sites"

# empties directory first
find "representative-binding-sites" -type f -delete

while IFS= read -r repr; do
	mid=${repr:4:2}

	cp "filtered-binding-sites/$mid/$repr" "representative-binding-sites/$repr"

done < cluster_reprs.txt
