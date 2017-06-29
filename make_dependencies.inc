alphadisk.o     : cgs.o globals.o grid.o rk4integr.o 
approximation.o : cgs.o globals.o 
balance.o       : threshold.o globals.o 
globals.o       : abundance.o cgs.o 
model_m1.o      : balance.o findzer.o alphadisk.o settings.o eulerintegr.o     \
    globals.o rk4settings.o threshold.o rk4integr.o fileunits.o cgs.o 
model_ss73.o    : eulerintegr.o globals.o grid.o rk4settings.o threshold.o     \
    rk4integr.o fileunits.o cgs.o 
diskvert-radial-1.o: model_ss73.o globals.o approximation.o settings.o cgs.o   \
    fileunits.o 
diskvert-m1.o   : globals.o grid.o settings.o model_m1.o fileunits.o results.o 
diskvert-ss73.o : model_ss73.o fileunits.o globals.o results.o settings.o 
query.o         : globals.o c_f_string.o 
relaxation.o    : fileunits.o globals.o cgs.o 
relaxation_C_interfaces.o: relaxation.o 
settings.o      : globals.o 
