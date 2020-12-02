'''
This file is used to launch the SiPM simulation using multiprocessing
module in order to speed up the process. It is possible to save ROOT
files of the futures extracted from the signals and the waveforms themselves.

Since the signals are kept in RAM untill they are saved on disk it is
recomended to launch this script on a small dataset if using the -W option
to save the waveforms. If you need to save waveforms of big datasets use the
file called 'wavedump.py'.
'''
from main import *

# Openig file
fname = args.txtfile


# Reading txt file
tot = 1
INPUT = []
OUTPUT = np.empty((0, 5))
OTHER = np.empty((0, 3))
with open(fname) as f:
    print(f'Opening file: {fname}')
    for line in f:
        line = line.strip().split(' ')
        eid = int(line[0])
        fid = int(line[1])
        fty = int(line[2])
        times = list(map(float, line[3:]))
        other = [eid, fid, fty]
        INPUT.append((times, other))
        if tot % 2000000 == 0:
            NFIB = len(INPUT)
            NEVT = eid
            print('\n===> Starting simulation <===\n')

            pool = Pool(args.jobs, initializer=initializeRandomPool, maxtasksperchild=10)
            res = pool.starmap_async(SiPM, INPUT, chunksize=100000)
            pool.close()
            pool.join()

            # Setting up results arrays
            output = np.empty(shape=(NFIB, 5), dtype='float32')
            other = np.empty(shape=(NFIB, 3), dtype='float32')

            for i, r in enumerate(res.get()):
                output[i, :] = r[0]
                other[i, :] = r[1]
            OUTPUT = np.vstack((OUTPUT, output))
            OTHER = np.vstack((OTHER, other))
            INPUT = []
        tot += 1

    NFIB = len(INPUT)
    NEVT = eid
    print('\n===> Starting simulation <===\n')

    pool = Pool(args.jobs, initializer=initializeRandomPool, maxtasksperchild=10)
    res = pool.starmap_async(SiPM, INPUT, chunksize=100000)
    pool.close()
    pool.join()

    # Setting up results arrays
    output = np.empty(shape=(NFIB, 5), dtype='float32')
    other = np.empty(shape=(NFIB, 3), dtype='float32')

    for i, r in enumerate(res.get()):
        output[i, :] = r[0]
        other[i, :] = r[1]
    OUTPUT = np.vstack((OUTPUT, output))
    OTHER = np.vstack((OTHER, other))


print('\n===> Simulation finished <===\n')

if args.graphics:
    print('\n===> Generating plots <===\n')
    somestats(OUTPUT)

np.savez_compressed(fname[:-5], output=OUTPUT, other=OTHER)
