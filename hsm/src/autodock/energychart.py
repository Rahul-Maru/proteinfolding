import csv
import json
from matplotlib import pyplot as plt
import numpy as np


def mid_avg(dat):
   
    q1 = np.percentile(dat, 25)
    q3 = np.percentile(dat, 75)
    
    return np.average([v for v in dat if q1 <= v <= q3])

# simple gaussian KDE
def kde_curve(data, num_points=200):
    data = np.asarray(data)
    n = len(data)
    
    std = np.std(data)
    if std == 0:
        std = 1e-6
    bw = 1.06 * std * (n ** (-1/5))
    xmin, xmax = np.min(data), np.max(data)
    if xmin == xmax:
        xmin -= 1
        xmax += 1
    x = np.linspace(xmin, xmax, num_points)
    # Gaussian kernel density estimate
    diffs = (x[None, :] - data[:, None]) / bw
    kernel_vals = np.exp(-0.5 * diffs**2) / (np.sqrt(2 * np.pi) * bw)
    y = kernel_vals.mean(axis=0)
    return x, y


with open("hsm/outs/autodock/energies.csv") as f:
	reader = csv.reader(f)
	lines = [row for row in reader]
	rows = {row[0]: float(row[1].strip()) for row in lines if row[1].strip()}
	no_energy = [row[0] for row in lines if not row[1].strip()]

prots = {k: v for k,v in sorted(rows.items(), key=lambda x : x[1]) if 'HSM' not in k}
energies = np.array(list(prots.values()))

hsm = {k: v for k,v in sorted(rows.items(), key=lambda x : x[1]) if 'HSM' in k} 
h_energies = np.array(list(hsm.values()))

print("# total: ", len(prots) + len(hsm))
print("# non-hsm: ", len(prots))
print("# hsm  ", len(hsm))
print("non-hsm avg energy: ", np.average(energies))
print("non-hsm avg energy (middle 50%): ", mid_avg(energies))
print("hsm avg energy: ", np.average(h_energies))
print("hsm avg energy (middle 50%): ", mid_avg(h_energies))

cutoff = -6
print(f"number of sites with energy < {cutoff}: ", len([v for v in prots.values() if v < cutoff]))

json.dump([k for k, v in prots.items() if v < cutoff], open("hsm/outs/autodock/prots_low_energy.json", "w"))

all_energies = np.concatenate([abs(energies), abs(h_energies)])
min_x = min(all_energies) - 0.3
max_x = max(all_energies) + 0.3

plt.subplot(1, 2, 1)
plt.hist(abs(energies), 20, color="xkcd:ultramarine blue", alpha=0.45, density=True)
x1, y1 = kde_curve(np.abs(energies))
plt.plot(x1, y1, color="xkcd:ultramarine blue", linewidth=2)

plt.xlim(min_x, max_x)
plt.xlabel("|E|")
plt.ylabel("Density (candidate sites)")

plt.subplot(1, 2, 2)
plt.hist(abs(h_energies), 5, color="xkcd:amber", alpha=0.45, density=True)
x2, y2 = kde_curve(np.abs(h_energies))

plt.xlim(min_x, max_x)
plt.plot(x2, y2, color="xkcd:amber", linewidth=2)
plt.xlabel("|E|")
plt.ylabel("Density (histamine sites)")

plt.show()