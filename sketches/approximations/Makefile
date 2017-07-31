
DATFILES = MCFK.dat MWFK.dat MQFK.dat MCFX.dat MWFX.dat MDFX.dat
COOLPNG = $(DATFILES:%.dat=%.cool2D.png)
DVFLAGS = -no-bf

all: $(DATFILES) $(COOLPNG) porown.png

%.dat: input.par

MCFX.dat:
	cat input.par | dv-mag-rx $(DVFLAGS) -corona -o MCFX
MWFX.dat:
	cat input.par | dv-mag-rx $(DVFLAGS) -compton -o MWFX
MDFX.dat:
	cat input.par | dv-mag-rx $(DVFLAGS) -o MDFX
MCFK.dat:
	cat input.par | dv-mag $(DVFLAGS) -corona -o MCFK
MWFK.dat:
	cat input.par | dv-mag $(DVFLAGS) -compton -o MWFK
MQFK.dat:
	cat input.par | dv-mag $(DVFLAGS) -compton2 -o MQFK

%.cool2D.png: %.dat
	diskvert-cooling2D --dpi 144 $< -o $@

porown.png: $(DATFILES)
	python plot.py

clean:
	$(RM) *.dat *.txt *.col
