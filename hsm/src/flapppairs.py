import json


dr="hsm/tools/MAPP-3D/MultipleSiteAlignment/flappsites"
pairs = json.load(open("hsm/outs/FLAPP/alignedsites.txt"))

l = "\n".join([f"{k}\t{s}" for k,v in pairs.items() for s in v]) + "\n"
print(l)

with open("hsm/tools/MAPP-3D/MultipleSiteAlignment/FPairs.txt", 'w') as f:
	f.write(l)
