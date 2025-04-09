import os
import subprocess

def filter_bsites():
	bsites = os.listdir("hsm/bsites")
	bsites.sort()
	bchains = os.listdir("hsm/bchains_final")
	bchains.sort()
	print(bsites)
	i = 0
	for s in bsites:
		if s not in bchains:
			print(s)
			i += 1
			subprocess.run(['rm', f'hsm/bsites_final/{s}'])
	print(i)


if __name__ == "__main__":
	filter_bsites()