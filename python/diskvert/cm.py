
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt

def mygrad(s,yl,yh,xg,yg,wg,xl,xh,sm):
    return LinearSegmentedColormap(s, {
        'red': (
            (0.0,yl,yl),
            (xl,0.0,0.0),
            (0.5,1.0,1.0),
            (1-xh,1.0-sm,1.0-sm),
            (1.0,yh,yh),
        ),
        'green': (
            (0.0,yl,yl),
            (0.5-xg,yg,yg),
            (0.5-0.5*xg, yg+wg*(1-yg), yg+wg*(1-yg)),
            (0.5,1,1),
            (0.5+0.5*xg, yg+wg*(1-yg), yg+wg*(1-yg)),
            (0.5+xg,yg,yg),
            (1.0,yl,yl),
        ),
        'blue': (
            (0.0,yh,yh),
            (xh,1.0-sm,1.0-sm),
            (0.5,0.97,0.97),
            (1-xl,0.0,0.0),
            (1.0,yl,yl),
        ),
    }, N = 1024)

plt.register_cmap(cmap=mygrad('bil1',0.1,0.36, 0.27,0.09,0.55, 0.4,0.25,0.04))
plt.register_cmap(cmap=mygrad('bil2',0.18,0.32, 0.3,0.25,0.75, 0.2,0.45,0.0))
