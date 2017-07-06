# Dominik Gronkiewicz (c) 2017
# this source is distributed under MIT License
# gronki@camk.edu.pl

import re

class pyminiconf(object):
    def __init__(self,f):
        buf = f.read() + "\n"

        # multiline strings
        rp = r'([a-zA-Z0-9\_\-\:]+)(?:[ ]+|[ ]*\=[ ]*)\"([^\"]*)\"\s+'
        for k,v in re.findall(rp,buf):
            dat[k] = v
        buf = re.sub(rp,'',buf)

        # komentarze
        buf = re.sub(r'\#[^\n]+','',buf)

        for k,v in re.findall(r'([a-zA-Z0-9\_\-\:]+)(?:[ ]+|[ ]*\=[ ]*)([\+\-]?[0-9]+)\s+',buf):
            setattr(self,k,int(v))
        for k,v in re.findall(r'([a-zA-Z0-9\_\-\:]+)(?:[ ]+|[ ]*\=[ ]*)([\+\-]?[0-9]+\.[0-9]+)\s+',buf):
            setattr(self,k,float(v))
        for k,v in re.findall(r'([a-zA-Z0-9\_\-\:]+)(?:[ ]+|[ ]*\=[ ]*)([\+\-]?[0-9]+\.?[0-9]*[eE][\+\-]?[0-9]+)\s+',buf):
            setattr(self,k,float(v))

        for k,v in re.findall(r'([a-zA-Z0-9\_\-\:]+)(?:[ ]+|[ ]*\=[ ]*)([yYtT]|[fFnN])\s+',buf):
            setattr(self,k,(v in ['T','t','Y','y']))

        for k,v in re.findall(r'([a-zA-Z0-9\_\-\:]+)(?:[ ]+|[ ]*\=[ ]*)([^0-9\-\+\s][^\s\#]+)\s+',buf):
            setattr(self,k,v)
