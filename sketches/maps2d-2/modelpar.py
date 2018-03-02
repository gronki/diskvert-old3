
class ModelDataset:
    def __init__(self, data):
        self._data = data
        self._okay = lambda d,k: d != None and hasattr(d,k) and type(getattr(d,k)) != str and d.converged
    def __getitem__(self,k):
        from numpy import vectorize
        return vectorize(lambda d: getattr(d,k) if self._okay(d,k) else np.nan)(self._data)
    def masked(self,k):
        from numpy import vectorize
        from numpy.ma import masked_where
        dm = vectorize(lambda d: not self._okay(d,k))(self._data)
        return masked_where(dm, self.__getitem__(k))

class ModelGrid:
    def __init__(self, mbh, mdot, radius, alpha, eta = 0.5, nu = 0.0):
        def fix(v):
            if hasattr(v, '__len__'):
                return list(v)
            else:
                return [v]
        self._keys = [ 'mbh', 'mdot', 'radius', 'alpha', 'eta', 'nu' ]
        self._vals = [ fix(v) for v in [ mbh, mdot, radius, alpha, eta, nu ] ]
        print self._vals
        self._lengths = [ len(x) for x in self._vals ]
        self._data_shape = filter(lambda n: n > 1, self._lengths)
        ix = range(len(self._keys))
        self._iter_fields = [ i for i,l in zip(ix,self._lengths) if l > 1 ]
        self._iter_axes = [ self._iter_fields.index(i) if i in self._iter_fields else None for i in ix ]

        self._flabels = [ 'H', 'M', 'R', 'A', 'B', 'N' ]
        self._llabels = [ '$M_{{\\rm BH}}$', '$\\dot{{m}}$',
            '$R / R_{{\\rm schw}}$', '$\\alpha_{{\\rm B}}$', '$\\eta$',
            '$\\nu_{{\\rm rec}}$' ]

        self._fn = lambda *ix: ''.join([ '{}{:02d}'.format(self._flabels[j1], j2 + 1) for j1,j2 in zip(self._iter_fields, ix) ])

    def axis(self,k):
        from numpy import array
        i = self._keys.index(k)
        return self._llabels[i], array(self._vals[i])

    def array_par(self):
        ix = range(len(self._keys))
        from numpy import ndarray, nditer
        data = ndarray(self._data_shape,
            dtype = [(l,'f4') for l in self._keys])
        for i,k,ia in zip(ix, self._keys, self._iter_axes):
            it = nditer(data, flags = ['multi_index'], op_flags = ['readwrite'])
            while not it.finished:
                it[0][k] = self._vals[i][it.multi_index[ia]] if ia != None \
                    else self._vals[i][0]
                print self._fn(*it.multi_index)
                it.iternext()
        return data

    def __len__(self):
        L = 1
        for v in self._vals: L *= len(v)
        return L

    def files(self):
        from itertools import product
        l = [range(i) for i in self._data_shape]
        for ix in product(*l):
            yield ix, self._fn(*ix)

    def parfiles(self):
        for ix,fn in self.files():
            par = {}
            for i,k,ia in zip(range(len(self._keys)), self._keys, self._iter_axes):
                par[k] = self._vals[i][ix[ia]] if ia != None \
                    else self._vals[i][0]
            fc = u"mbh    {mbh:8.4e}\n"    \
                  "mdot   {mdot:8.6f}\n"   \
                  "alpha  {alpha:8.6f}\n"  \
                  "radius {radius:8.3f}\n" \
                  "zeta   {eta:8.6f}\n"    \
                  "nu     {nu:8.3f}\n".format(**par)
            yield ix,fn,fc

    def write_par(self, path = lambda fn: 'data/{}.par'.format(fn)):
        for ix,fn,fc in self.parfiles():
            with open(path(fn),'w') as f: f.write(fc)

    def read_dataset(self, path = lambda fn: 'data/{}.txt'.format(fn)):
        from diskvert import pyminiconf
        from numpy import ndarray
        data = ndarray(self._data_shape, dtype = object)
        for ix, fn in self.files():
            try:
                with open(path(fn), 'r') as f: data[tuple(ix)] = pyminiconf(f)
            except IOError:
                data[tuple(ix)] = None
        return ModelDataset(data)

    # TODO: axes2d lub cos takeico co zwraca 2 osie

if __name__ == '__main__':
    import numpy as np
    gd = ModelGrid(mbh = 10, mdot = [0.01,0.03,0.1,0.3], alpha = [0.01,0.03,0.1], radius = 10)
    # gd.write_par(lambda fn: 'data/{}.par'.format(fn))
    D = gd.read_dataset(lambda fn: 'data/{}Q.txt'.format(fn))
    print D['tavg_tmin']
    # xl, xv = gd.axis('mdot')
    # yl, yv = gd.axis('alpha')
    # plt.pcolor(xv, yv, D['tcor'])
