alphadisk.o     : rk4integr.o eulerintegr.o ss73solution.o globals.o grid.o    \
    rk4settings.o threshold.o heatbil.o fileunits.o cgs.o 
alphasimp.o     : cgs.o globals.o grid.o rk4integr.o fileunits.o 
balance.o       : threshold.o globals.o 
globals.o       : cgs.o fileunits.o 
heatbil.o       : cgs.o globals.o 
modelmag.o      : balance.o findzer.o settings.o eulerintegr.o alphasimp.o     \
    globals.o rk4settings.o threshold.o rk4integr.o fileunits.o cgs.o 
dv-alpha.o      : fileunits.o alphadisk.o globals.o results.o settings.o 
dv-alpha-rx.o   : ss73solution.o globals.o grid.o settings.o fileunits.o       \
    relaxation.o 
dv-mag.o        : globals.o grid.o settings.o fileunits.o modelmag.o results.o 
dv-mag-rx.o     : ss73solution.o globals.o grid.o settings.o fileunits.o       \
    relaxation.o cgs.o 
dv-rad1.o       : ss73solution.o globals.o settings.o alphadisk.o cgs.o        \
    fileunits.o 
query.o         : globals.o c_f_string.o 
relaxation.o    : fileunits.o globals.o cgs.o 
relaxation_c.o  : relaxation.o 
settings.o      : globals.o 
ss73solution.o  : cgs.o globals.o 
