from collections import defaultdict
import json

cutoff = float(open("hsm/src/autodock/cutoff").read().strip())
clust_cutoff = 10

prots = json.load(open("hsm/outs/autodock/prots_low_energy.json"))

clusts = defaultdict(list)

for p in prots:
    with open(f"hsm/outs/autodock/results/{p}/1.dlg") as f:
        for l in f:
            # iterate through clustering histogram
            if l.startswith("_____|___________|_____|___________|_____|____:____|____:____|____:____|____:___"):
                for l2 in f:
                    if l2.startswith("_____|___________|_____|___________|_____|______________________________________"):
                        break
    
                    # 0: cluster rank, 1: min binding energy, 2: run no,
                    #   3: mean binding energy, 4: number in cluster
                    toks = [float(k) for k in l2.split('|')[:-1]]

                    # if mean binding energy of rank one cluster ≥ cutoff, reject the site
                    if toks[0] == 1 and toks[3] >= cutoff:
                        break

                    # store the cluster data
                    clusts[p].append(toks)

                else:
                    raise RuntimeError(f"histogram end separator not found for {p}")

                break

        # print()

final_clusts = defaultdict(list)
for p, c in clusts.items():
    print(p, c)
    for r in range(len(c)):
        # if mean binding energy of rank r >= cutoff, reject the site
        if c[r][3] >= cutoff:
            break
        # if number in cluster >= clust_cutoff, add the cluster to the final clusters
        if c[r][4] >= clust_cutoff:
            # rank and run #
            final_clusts[p].append((int(c[r][0]), int(c[r][2])))

print('\n'.join([f"{k}: {v}" for k, v in final_clusts.items()]))            
print(sum([len(v) for _, v  in final_clusts.items()]))
print(len(final_clusts))

json.dump(final_clusts, open("hsm/outs/autodock/docking_clusts.json", "w"))

# print('\n'.join([f"{k}: {v}" for k, v in nums.items()]))
# print(len(clusts))