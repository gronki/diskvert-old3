#!/usr/bin/env python
# coding: utf-8

files = {}
commands = []
commands_pl2d = []


cml1 = dict(
    AK = 'dv-alpha',
    AX = 'dv-alpha-rx -n 300',
    MK = 'dv-mag -top 100',
    MX = 'dv-mag-rx -top 37 -n 300',
)
cml2 = dict(
    D = '',
    C = '-corona',
    W = '-compton',
)
cml3 = dict(
    F = '-no-bf',
    FB = '',
)

for ss in ['X','K']:
    for sm in ['D','W','C']:
        for sh in ['A','M']:
            for so in ['F','FB']:
                if ss == 'X' and sm == 'W': continue
                if sh == 'A' and sm != 'D': continue
                name = sh + sm + so + ss
                cmdline = " ".join([cml1[sh+ss],cml2[sm],cml3[so],'-o',name])
                commands.append("cat input.par | {0} &".format(cmdline))
                if ss == 'K' and sh == 'M':
                    commands_pl2d.append("tar czf {0}.tgz {0}.{{dat,col,txt}}".format(name))
                    commands_pl2d.append("diskvert-cooling2D {0}.tgz".format(name))
                files[name] = sh+ss
        commands.append("wait")

if __name__ == '__main__':
    print "#!/usr/bin/env bash\n\nrm -f *.{png,col,dat,txt,log,tgz}\n\n"
    for c in commands: print c
    print 'python plot.py'
    for c in commands_pl2d: print c
