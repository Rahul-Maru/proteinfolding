import csv
from matplotlib import pyplot as plt
import ast

with open("hsm/reprgrid.csv") as f:
	r = csv.reader(f, delimiter=';')
	lines = [[k.strip() for k in l] for l in r][:-1]
print (lines)

data_by_class = {}
for line in lines:
	try:
		x, y, cls = float(line[0]), float(line[1]), ast.literal_eval(line[2])[0]
	except Exception as e:
		print(line)
		raise e
	if cls not in data_by_class:
		data_by_class[cls] = {'x': [], 'y': []}
	data_by_class[cls]['x'].append(x)
	data_by_class[cls]['y'].append(y)

plt.figure(figsize=(10, 6))
for cls, data in data_by_class.items():

	scatter = plt.scatter([], [])  # Create temporary scatter to get color
	color = scatter.get_facecolors()[0]
	scatter.remove()
	
	# Plot 1x1 squares centered at each point with matching color
	for x, y in zip(data['x'], data['y']):
		square = plt.Rectangle((x - 0.025, y - 0.5), 0.05, 1, alpha=0.42, color=color)
		plt.gca().add_patch(square)

	plt.scatter(data['x'], data['y'], label=cls, s=50, color=color)

plt.xlabel('M_dist min score cutoff')
plt.ylabel('Min # of common residues')
plt.title('Classification Scatter Plot')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()