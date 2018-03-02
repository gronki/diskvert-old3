from sys import argv
from math import sqrt

data = []
for fn in argv[1:]:
    with open(fn) as f:
        data.append(float(f.read()))

n = len(data)
mean = sum(data) / n
stdev = sqrt(sum([ (d - mean)**2 for d in data ]) / (n**2 - n))

print u'{:8.3f} {:8.3f}'.format(mean,stdev)
