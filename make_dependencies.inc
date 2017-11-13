alphadisk.o     : rk4integr.o eulerintegr.o ss73solution.o globals.o grid.o    \
    rk4settings.o threshold.o heatbil.o fileunits.o cgs.o 
balance.o       : threshold.o globals.o 
globals.o       : cgs.o 
heatbil.o       : cgs.o globals.o 
modelmag.o      : balance.o findzer.o settings.o eulerintegr.o ss73solution.o  \
    globals.o rk4settings.o threshold.o heatbil.o fileunits.o cgs.o            \
    rk4integr.o 
dv-alpha.o      : fileunits.o alphadisk.o globals.o results.o settings.o 
dv-alpha-rx.o   : settings.o relaxutils.o rxsettings.o ss73solution.o          \
    globals.o grid.o relaxation.o fileunits.o 
dv-mag.o        : globals.o grid.o settings.o fileunits.o modelmag.o results.o 
dv-mag-rx.o     : settings.o relaxutils.o rxsettings.o ss73solution.o          \
    globals.o grid.o relaxation.o interpol.o fileunits.o deriv.o 
dv-rad1.o       : ss73solution.o globals.o settings.o alphadisk.o cgs.o        \
    fileunits.o 
relaxation.o    : fileunits.o globals.o cgs.o 
relaxation_c.o  : relaxation.o 
relaxutils.o    : globals.o 
rk4settings.o   : globals.o 
settings.o      : fileunits.o globals.o 
ss73solution.o  : cgs.o globals.o 
