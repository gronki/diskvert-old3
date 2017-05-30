# Dominik Gronkiewicz  2016
# gronki@camk.edu.pl

import re

def pyminiconf_read(f):
    buf = f.read() + "\n"
    dat = {}

    # multiline strings
    rp = r'([a-zA-Z0-9\_\-\:]+)([ ]+|[ ]*\=[ ]*)\"([^\"]*)\"\s+'
    for k,s,v in re.findall(rp,buf):
        dat[k] = v
    buf = re.sub(rp,'',buf)

    # komentarze
    buf = re.sub(r'\#[^\n]+','',buf)

    for k,s,v in re.findall(r'([a-zA-Z0-9\_\-\:]+)([ ]+|[ ]*\=[ ]*)([\+\-]?[0-9]+)\s+',buf):
        dat[k] = int(v)
    for k,s,v in re.findall(r'([a-zA-Z0-9\_\-\:]+)([ ]+|[ ]*\=[ ]*)([\+\-]?[0-9]+\.[0-9]+)\s+',buf):
        dat[k] = float(v)
    for k,s,v in re.findall(r'([a-zA-Z0-9\_\-\:]+)([ ]+|[ ]*\=[ ]*)([\+\-]?[0-9]+\.?[0-9]*[eE][\+\-]?[0-9]+)\s+',buf):
        dat[k] = float(v)

    for k,s,v in re.findall(r'([a-zA-Z0-9\_\-\:]+)([ ]+|[ ]*\=[ ]*)([^0-9\-\+\s][^\s\#]+)\s+',buf):
        dat[k] = v

    return dat
