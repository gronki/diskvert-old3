# Estimation of model truncation height

This mini project is to estimate the upper boundary of the computational interval. For highly magnetized models, this height must be large, while for the models when reconnection parameter is significant, extending the range to extremely low densities will result in numerical instabilities

## Tutorial

```bash
mkdir data
python par.py
find data -name \*.par | parallel --eta cat {} \| dv-mag-rx -o {.}
find data -name \*.dat | parallel python plot2.py
```
