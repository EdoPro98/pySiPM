from sipm import *
import time

NPE = 10
NEVTS = 50000

pool = multiprocessing.Pool(args.jobs, initializer=initializeRandomPool, initargs=(0,))

times = []
other = []

for i in range(NEVTS):
    n = np.random.poisson(NPE)
    times.append(np.abs(np.random.normal(20, 5, n)).tolist())
    other.append([])

input = zip(times, other)c
startingtime = time.time_ns()
results = pool.starmap_async(SiPM, input)
pool.close()
pool.join()
endingtime = time.time_ns()

elapsedtime = (endingtime - startingtime)
eventsperms = NEVTS / (elapsedtime / 1e6)
print(f'\nElapsed time: {elapsedtime / 1e6:.2f} ms [{elapsedtime / 1e9:.2f} s]')
print(f'Simulation speed: {eventsperms:.3f} signals/ms [{1/eventsperms:.2f} ms per signal]')

features = np.empty((NEVTS, 5), dtype='float')
eventinfo = []
for i, r in enumerate(results.get()):
    features[i, :] = r[0]
    eventinfo.append(r[1])
    signal = r[2]


if args.graphics:
    somestats(features)
