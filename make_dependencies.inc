alphadisk.o     : cgs.o globals.o precision.o 
balance.o       : threshold.o globals.o precision.o 
globals.o       : kramers.o rk4integr.o cgs.o precision.o space.o 
model_m1.o      : balance.o findzer.o alphadisk.o settings.o precision.o       \
    eulerintegr.o kramers.o globals.o threshold.o cgs.o rk4integr.o 
model_ss73.o    : globals.o threshold.o rk4integr.o cgs.o space.o precision.o  \
    eulerintegr.o 
diskvertp.o     : threshold.o cgs.o setup.o globals.o precision.o 
disk-properties.o: setup.o globals.o precision.o settings.o 
diskvert-m1.o   : globals.o settings.o model_m1.o setup.o precision.o          \
    results.o 
diskvert-ss73.o : model_ss73.o globals.o settings.o setup.o precision.o        \
    results.o 
coefficients.o  : cgs.o globals.o precision.o 
relaxation.o    : lapack_interfaces.o cgs.o globals.o precision.o              \
    coefficients.o 
relaxation_C_interfaces.o: precision.o relaxation.o 
results.o       : precision.o 
settings.o      : globals.o results.o 
setup.o         : globals.o settings.o 
bisect.o        : precision.o 
cgs.o           : precision.o 
deriv.o         : precision.o 
eulerintegr.o   : precision.o 
findzer.o       : precision.o 
findzer_multi.o : precision.o 
histogram.o     : precision.o 
interpol.o      : precision.o 
kramers.o       : cgs.o precision.o 
linsect.o       : precision.o 
random.o        : precision.o 
rk4integr.o     : precision.o 
space.o         : precision.o 
threshold.o     : precision.o 
