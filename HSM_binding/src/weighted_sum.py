import csv

def main():
	with open("HSM/MAPP-3D/MultipleSiteAlignment/mdist.csv") as f:
		reader = csv.reader(f)
		next(reader)

		prots = {}

		for row in reader:
			if row[0] not in prots:
				prots[row[0]] = 0
			if row[1] not in prots:
				prots[row[1]] = 0

			prots[row[0]] += float(row[2])
			prots[row[1]] += float(row[2])
		
		print(prots)
		print(max(prots, key=prots.get))
			

if __name__ == "__main__":
	main()