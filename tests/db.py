import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import sys
from scipy.optimize import curve_fit as fit
from uncertainties import ufloat
import numexpr as ne
import math
sys.path.append('../')
from sipm import *

N = 5000

PEAKS = np.empty(N, dtype=np.float32)
SIGMA = []
EXPECTED_DCR = DCR * INTGATE
EXPECTED_XT = DCR * INTGATE * XT

args.debug = True

for i in range(N):
    out, other, signal, debug = SiPM([20]*np.random.randint(0, 2), ())
    PEAKS[i] = out[0]
    if (debug[0] + debug[1]) == 0:
        SIGMA.append(np.std(signal))
SIGMA = np.asarray(SIGMA)

h = np.histogram(PEAKS, 300, density=True, range=(0, 2))
x = h[1]
x = (x[1:] + x[:-1])/2
y = h[0]

x_one_peak = x[(x > 0.5) & (x < 1.5)]
y_one_peak = y[(x > 0.5) & (x < 1.5)]


def gaus(x, a, mu, sigma):
    return ne.evaluate("a*exp(-0.5*((mu - x) / sigma)**2)")


par1, cov1 = fit(gaus, x_one_peak, y_one_peak, bounds=(0, np.inf))
err1 = np.diag(cov1)

signal = ufloat(par1[1], par1[2])
noise = ufloat(np.mean(SIGMA), np.std(SIGMA)/SIGMA.size**0.5)
snr = 20 * math.log10(par1[1] / np.mean(SIGMA))
print(f"Snr is {snr:.1f} dB")

plt.plot(x_one_peak, gaus(x_one_peak, *par1), 'r')
hep.histplot(h, histtype='fill', color='k')
plt.show()
