from sys import argv, stderr, stdout
import numpy as np

data = []
for fn in argv[1:]:
    with open(fn) as f:
        data.append(float(f.read()))

times = np.ma.masked_array(data)
if times.count() >= 5:
    times[times.argmin()] = np.ma.masked
    times[times.argmax()] = np.ma.masked

stderr.write(str(data) + "\n")

n = times.count()
mean = np.sum(times) / n
stdev = np.sqrt(np.sum((times - mean)**2) / (n - 1))

print(u'{:8.3f} {:8.3f}\n'.format(mean,stdev))
