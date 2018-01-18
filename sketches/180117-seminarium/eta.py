from par import *
from matplotlib.cm import coolwarm as cm
from matplotlib.cm import coolwarm_r as cm_r
import matplotlib.pyplot as plt

betas = np.logspace(-1, 3, 80)
etas = [0.01, 0.03, 0.1, 0.31, 1.0]

plt.figure()
for i in range(len(etas)):
    alphas = 2 * etas[i] / (1 + betas)
    plt.plot(betas, alphas, label = '$\\eta$ = {:.2f}'.format(etas[i]), color = cm(i / 4.0))

plt.plot(betas, 0.25 / np.sqrt(betas), '--', color = '#AA9AC9', label = 'modeling')

plt.loglog()
plt.legend()
plt.xlabel('$\\beta_0$')
plt.ylabel('$\\alpha_B$')
plt.tight_layout()
plt.savefig('etas.png')
