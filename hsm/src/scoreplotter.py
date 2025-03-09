import csv
from matplotlib import pyplot as plt
import numpy as np

def main():
	with open("hsm/outs/Clustal/network.csv") as f:
		r = csv.reader(f)
		next(r)
		dat = [i[2] for i in r]
	
	n = len(dat)
	x = np.linspace(0, 1, 101)
	y = [len([j for s in dat if (j:=float(s))<=i])/n for i in x]

	plt.plot(x, y, label="cumulative frequency of scores")

	plt.show()

	

if __name__ == "__main__":
	main()