import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import sys
from uncertainties import ufloat
sys.path.append('../')
from sipm import *

N = 100000
TRUEDCR = np.empty(N,dtype=np.uint8)
TRUEXT = np.empty(N,dtype=np.uint8)
PEAKS = np.empty(N,dtype=np.float32)
INTEGRALS = np.empty(N,dtype=np.float32)
EXPECTED_DCR = DCR * INTGATE
EXPECTED_XT = DCR * INTGATE * XT

args.NAP = True
args.debug = True
times = [[]]*N

pool = Pool(4, initializer=initializeRandomPool)
mapresult = pool.map_async(SiPM,times)
pool.close()
pool.join()

for i,result in enumerate(mapresult.get()):
    out, signal, other, debug = result
    PEAKS[i] = out[0]
    INTEGRALS[i] = out[1]
    TRUEDCR[i] = debug[1]
    TRUEXT[i] = debug[2]

pedestal = np.mean(PEAKS[PEAKS < 0.5])
MEASUREDDCR = np.uint8(PEAKS > 0.5 + pedestal)
MEASUREDXT = np.uint8((PEAKS > 1.5 + pedestal) & (PEAKS < 2.5 + pedestal))

mctruth_mean_dcr = ufloat(TRUEDCR.sum(), TRUEDCR.sum()**0.5)
mctruth_mean_xt = ufloat(TRUEXT.sum(), TRUEXT.sum()**0.5)

measured_mean_dcr = ufloat(MEASUREDDCR.sum(), MEASUREDDCR.sum()**0.5)
measured_mean_xt = ufloat(MEASUREDXT.sum(), MEASUREDXT.sum()**0.5)

mctruth_xt = mctruth_mean_xt / mctruth_mean_dcr
mctruth_dcr = mctruth_mean_dcr / SIGLEN * 1e9 / N
measured_xt = measured_mean_xt / measured_mean_dcr
measured_dcr = measured_mean_dcr / (INTGATE * SAMPLING) * 1e9 / N

print(f"Expected DCR: {DCR/1e3:.3f} kHz")
print(f"MC Truth DCR {mctruth_dcr/1e3:.3f} kHz")
print(f"Measured DCR {measured_dcr/1e3:.3f} kHz")
print()
print(f"Expected XT: {XT*100:.3f} %")
print(f"MC Truth XT {mctruth_xt*100:.3f} %")
print(f"Measured XT {measured_xt*100:.3f} %")

fig, ax = plt.subplots(2, 1)
h_peak = np.histogram(PEAKS, 500)
hep.histplot(h_peak, color='k', histtype='fill', ax=ax[0])
ax[0].set_yscale('log')

x = np.sort(PEAKS)
y = np.empty_like(x)
for i, p in enumerate(x):
    n = np.count_nonzero(PEAKS > p) / PEAKS.size
    rate = n / (INTGATE * SAMPLING) * 1e9
    y[i] = rate

ax[1].plot(x, y/1e3, 'k', label='Events over threshold')
ax[1].set_yscale('log')
ax[1].hlines(DCR/1e3,x[0],x[-1],'r', label='Expected DCR')
ax[1].hlines(XT * DCR/1e3,x[0],x[-1],'r', label='Expected XT')
plt.show()
