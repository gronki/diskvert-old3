from re import match

processor = -1
info = {}

with open('/proc/cpuinfo','r') as f:
    for l in f:
        m = match(r'^(processor|cpu family|model|model name|cpu MHz|cache size|core id|physical id)\s+:\s+([^\n]+)\n$', l)
        if not m: continue
        k,v = m.groups()
        if k == 'processor':
            processor = int(v)
            info[processor] = {}
        info[processor][k] = v
    for inf in info.values():
        print "{physical id}:{core id}:{processor} F{cpu family}M{model} {model name}\n    at {cpu MHz} MHz, cache {cache size}".format(**inf)
