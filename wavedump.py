'''
This file is used to launch the SiPM simulation and save the digitized waves of very large datasets.
Each SiPM event is simulated in a serial way and each wave is immediately
saved in a HDF5 datasets along with the features extracted and the
geometrical informations about the SiPM in the calorimeter.

launch with: python wavedump.py -W outfilename.hdf5

WARNING: make sure to have enough free space in the
working directory to save all the data!
'''
from sipm import *
from itertools import starmap

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

# Setting up list containing outputs and HDF5 files
fnameout = args.wavedump
hf = h5py.File(fnameout, 'w')
dset1 = hf.create_dataset('Waveforms', shape=(NFIB, SIGPTS), dtype='f', compression='gzip', compression_opts=9, chunks=(10, SIGPTS))
dset2 = hf.create_dataset('Features', shape=(NFIB, 5), dtype='f', compression='gzip', compression_opts=9)
dset3 = hf.create_dataset('Geometry', shape=(NFIB, 6), dtype='f', compression='gzip', compression_opts=9)
output = np.empty(shape=(NFIB, 5), dtype='float32')
other = np.empty(shape=(NFIB, 6), dtype='float32')


print('\n===> Starting simulation <===\n')
# Loop using multiprocessing on chunks of data
# Output is stored in HDF5 file inside the loop
pool = Pool(args.jobs, initializer=initializeRandomPool, maxtasksperchild=2**16)
CHUNK = 100000
if CHUNK > NFIB:
    CHUNK = NFIB
temp = 0
i = 0
Ts = time.time()
while i < NFIB:
    INPUT = zip(TIMES[i:i+CHUNK], OTHER[i:i+CHUNK])
    res = pool.starmap_async(SiPM, INPUT)
    res.wait()
    for j, r in enumerate(res.get()):
        dset1[i+j, :] = r[2]
        dset2[i+j, :] = r[0]
        dset3[i+j, :] = r[1]
        output[i+j, :] = r[0]
        other[i+j, :] = r[1]
    if OTHER[i][0] - temp > 10:
        hf.flush()
        print(f'Events simulated: {int(OTHER[i][0]):d} of {NEVT:d}')
        temp = OTHER[i][0]
    i += CHUNK
INPUT = zip(TIMES[i:], OTHER[i:])
res = starmap(SiPM, INPUT)
for j, r in enumerate(res):
    dset1[i+j, :] = r[2]
    dset2[i+j, :] = r[0]
    dset3[i+j, :] = r[1]
    output[i+j, :] = r[0]
    other[i+j, :] = r[1]
hf.close()
Te = time.time()

print('\n===> Simulation finished <===\n')
print(f'Execution time: {(Te-Ts):.2f}s')
print(f'Events per second: {NEVT/(Te-Ts):.2f}')

if args.write:
    print('\n===> Writing results <===\n')
    SaveFile(args.write, output, other)
