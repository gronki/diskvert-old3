alphadisk.o     : cgs.o globals.o grid.o 
diskapprox.o    : cgs.o globals.o 
balance.o       : threshold.o globals.o 
globals.o       : abundance.o cgs.o 
model_m1.o      : balance.o settings.o alphadisk.o findzer.o eulerintegr.o     \
    globals.o rk4settings.o threshold.o cgs.o rk4integr.o 
model_ss73.o    : globals.o grid.o threshold.o cgs.o rk4integr.o eulerintegr.o 
disk-properties.o: globals.o settings.o 
diskvert-radial-1.o: model_ss73.o cgs.o globals.o diskapprox.o settings.o 
diskvert-m1.o   : globals.o grid.o settings.o model_m1.o fileunits.o results.o 
diskvert-ss73.o : model_ss73.o globals.o results.o settings.o 
query.o         : globals.o c_f_string.o 
coefficients.o  : cgs.o globals.o 
relaxation.o    : cgs.o globals.o coefficients.o 
relaxation_C_interfaces.o: relaxation.o 
settings.o      : globals.o grid.o results.o 
