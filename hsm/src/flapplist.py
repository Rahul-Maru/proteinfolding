import os


nrsites = os.listdir("hsm/tools/FLAPP/sites")
hsmsites = [site for site in nrsites if 'HSM' in site]

# print(len(nrsites), len(hsmsites))
# print(len(nrsites)*len(hsmsites))

i = 0
for site in hsmsites:
	for site2 in nrsites:
		# i += 1
		# if not i % 50000:
		# 	print(i)
		print(f"{site}\t{site2}")

# reverse
# for site in hsmsites:
# 	for site2 in nrsites:
# 		print(f"{site}\t{site2}")
