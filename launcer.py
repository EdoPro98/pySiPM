'''
This file is used to launch the SiPM simulation using multiprocessing
module in order to speed up the process. It is possible to save ROOT
files of the futures extracted from the signals and the waveforms themselves.

Since the signals are kept in RAM untill they are saved on disk it is
recomended to launch this script on a small dataset if using the -W option
to save the waveforms. If you need to save waveforms of big datasets use the
file called 'wavedump.py'.
'''
from sipm import *

# Openig file
fname = 'Data/E40GeV_1kevents.txt'

# Reading txt file
with open(fname) as f:
    print(f'Opening file: {fname}')
    TIMES = []
    OTHER = []
    temp = 0
    for line in f:
        if not line.strip():
            continue
        L = line.split(' ')
        t = list(map(float,L[7:]))
        TIMES.append(t)
        OTHER.append((int(L[0]), int(L[1] == 'Scin'), int(L[2]), float(L[3]), float(L[4]), float(L[5])))
        if int(L[0]) % 100 == 0 and int(L[0]) > 0 and temp != L[0]:
            temp = L[0]
            print(f'Reading event: {int(L[0]):d}', end='\r')

OTHER = np.array(OTHER)
INPUT = zip(TIMES, OTHER)
NFIB = len(TIMES)
NEVT = int(L[0])
del TIMES
del OTHER

# Setting up results arrays
pool = Pool(processes=nJobs, initializer=initializeRandomPool,
            maxtasksperchild=16384)
output = np.empty(shape=(NFIB, 5), dtype='float32')
other = np.empty(shape=(NFIB, 6), dtype='float32')
if args.wavedump:
    signals = np.empty(shape=(NFIB, SIGPTS), dtype='float16')
print('\n===> Starting simulation <===\n')

# Launching simulation
Ts = time.time()
res = pool.starmap_async(SiPM, INPUT)
pool.close()
pool.join()
Te = time.time()

for i, r in enumerate(res.get()):
    output[i, :] = r[0]
    other[i, :] = r[1]
    if args.wavedump:
        signals[i, :] = r[2]

print('\n===> Simulation finished <===\n')
print(f'Execution time: {(Te-Ts):.2f}s')
print(f'Events per second: {NEVT/(Te-Ts):.2f}')

if args.graphics:
    print('\n===> Generating plots <===\n')
    somestats(output)

if args.write:
    print('\n===> Writing results <===\n')
    SaveFile(args.write, output, other)

if args.wavedump:
    print('\n===> Writing waveforms <===\n')
    SaveWaves(args.wavedump, signals)
