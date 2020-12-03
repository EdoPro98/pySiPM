from main import *

NPE = 10
NEVTS = 50000

pool = multiprocessing.Pool(args.jobs, initializer=initializeRandomPool, initargs=(0,))

times = []
other = []

for i in range(NEVTS):
    n = np.random.poisson(NPE)
    times.append(np.abs(np.random.normal(20, 5, n)).tolist())
    other.append([])

input = zip(times, other)
results = pool.starmap_async(SiPM, input)
pool.close()
pool.join()

features = np.empty((NEVTS, 5), dtype='float')
eventinfo = []
for i, r in enumerate(results.get()):
    features[i, :] = r[0]
    eventinfo.append(r[1])
    signal = r[2]


if args.graphics:
    somestats(features)
