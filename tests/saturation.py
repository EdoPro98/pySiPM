import sys
sys.path.append('../')
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as fit
import uncertainties as unc
import uncertainties.unumpy as unp
import numexpr as ne
from sipm import *

N = 1000
PEAKS = np.empty(N, dtype=np.float32)

print(f"Running saturation with {CELLSIZE:.0f} um SiPM")
TIMES = []
REALPE = []
for i in range(N):
    TIMES.append(np.random.normal(20, 0.01, np.random.randint(1, 2000)).tolist())
    REALPE.append(len(TIMES[-1]))
p = Pool(4)
mapresult = p.map_async(SiPM,TIMES)
p.close()
p.join()

REALPE = np.asarray(REALPE, dtype=np.uint16)

for i,result in enumerate(mapresult.get()):
    out, other, signal, debug = result
    PEAKS[i] = out[0]


def expo(x, a, b, c):
    return ne.evaluate("a * (1 - exp(-x/(a * b))) + c")


par, cov = fit(expo, REALPE, PEAKS, p0=[NCELL, 1, 0], bounds=(0, np.inf))
err = np.diag(cov)**0.5
print(f"Real number of cells: {NCELL+1}")
print(
    f"From fit: NCELL = {par[0]:.2f} +/- {err[0]:.2f}\tDPP = {par[1]:.2f} +/- {err[1]:.2f}")

low_x = REALPE[REALPE < 50]
low_y = PEAKS[REALPE < 50]
par_linear = np.polyfit(low_x, low_y, 1)

px = np.linspace(1,REALPE.max(),100)
plt.scatter(REALPE, PEAKS, c='k', s=3, label='Data')
plt.plot((0, REALPE.max()), np.polyval(par_linear, (0, REALPE.max())), 'b', label='Expected linearity')
plt.title(f"Response: pitch = {CELLSIZE:.0f} $\mu m$")
plt.plot(px, expo(px, *par), 'r', label='Occupancy fit')
plt.xlim(0, REALPE.max())
plt.ylim(0, REALPE.max())
plt.xlabel('Measured peak value')
plt.ylabel('Real number of photoelectrons')
plt.legend(fancybox=False, edgecolor='k')

fig, ax = plt.subplots(2,1,gridspec_kw={'height_ratios':[2,1]})
recope = -par[0] * np.log(1 - PEAKS * par[1] / par[0]) + par[2]
ax[0].scatter(REALPE, recope, c='k', s=3, label='Data')
ax[0].set_xlabel('Real number of photoelectrons')
ax[0].set_ylabel('Reconstructed number of photoelectrons')
ax[0].set_title(f"Reconstructed number of pe: pitch = {CELLSIZE:.0f} $\mu m$")
ax[0].plot((0, REALPE.max()), (0, REALPE.max()), 'b', label='Expected linearity')
plt.legend(fancybox=False, edgecolor='k')
ax[1].scatter(REALPE, 100*(REALPE - recope)/REALPE,c='k',s=2)
ax[1].fill_between((0,REALPE.max()),(-5,-5),(5,5),color='r',alpha=0.2,label="$\pm 5\%$")
ax[1].legend(fancybox=False, edgecolor='k')
plt.show()
