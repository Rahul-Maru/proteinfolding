from collections import defaultdict
import json


prots = json.load(open("hsm/outs/autodock/prots_low_energy.json"))

clusts = defaultdict(list)
nums = defaultdict(list)

for p in prots:
    with open(f"hsm/outs/autodock/{p}/1.dlg") as f:
        print(p, end=": \n")
        for l in f:
            if l.startswith("_____|___________|_____|___________|_____|____:____|____:____|____:____|____:___"):
                for l2 in f:
                    if l2.startswith("_____|___________|_____|___________|_____|______________________________________"):
                        break
    
                    print(l2, end="")
                    toks = [float(k) for k in l2.split('|')[:-1]]
                    if toks[0] == 1 and toks[3] > -6:
                        clusts.pop(p, None)
                        break

                    clusts[p].append(toks)
                    nums[p].append(toks[4])

                else:
                    raise

                break

        print()

print('\n'.join([f"{k}: {v}" for k, v in nums.items()]))
print(len(clusts))