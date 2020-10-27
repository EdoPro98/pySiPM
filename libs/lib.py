"""In this file I define all the functions I will use in the main file of simulation."""
# EDITING THIS FILE MAY SERIOUSLY COMPROMISE SIMULATION BEHAVIOUR!

from libs.FortranFunctions import rollfortran, signalgenfortran, sortfortran
from libs.FortranFunctions import frandom
from variables import *


def addDCR(rate):
    """!@brief Generation of dark count events."""
    """!
    Function that generates times for dark count events (DCR).
    Events are generated as an uniform poisson process, hence the time
    distance between consecutive DCR events follows an exponential
    distribution with a mean value of: @f$\tau=\frac{1}{DCR}@f$

    @param rate:        Dark counts rate un Hz.

    @return dcrTime:    List containing times of DCR events in ns.
    """
    dcrTime = []
    last = 0
    while last < SIGLEN:  # Generate from exponential distribution of delays
        last = frandom.randexp(1e9 / rate, 1)
        dcrTime.append(last.item())
    dcrTime.pop(-1)
    return dcrTime


def HitCells(n):
    """!@brief Generation of cell IDs."""
    """!
    Function that generated cell IDs for each photoelectron in the list.@n
    Cell IDs are integers in range [0 - NCELLS] and are generated randomly.
    An ID can appear multiple times, meaning that the corresponding cell has
    generated an avalanche due to a photoelectron multiple times.

    @param n:     Number of photoelectrons that have to be simulated.
    @return idx:    List containing the ID of each hitted cell.
    """
    idx = list(frandom.randint(NCELL, n))
    return idx


def addXT(times, idx, xt):
    """!@brief Generation of optical crosstalk events."""
    """!
    Function that generates optical crosstalk events(XT). Each photoelectron
    has generates a poissonian number of XT events accordingly to the crosstalk
    probability. XT events are generated in the 8 neighbouring cells and have
    the same time of the main cell.@n
    XT events can also generate other XT events in theyr neighbouring cells.

    @param times:   List containing the time of each SiPM cell avalanche,
                    including DCR events.
    @param idx:     List containing the ID of each fired cell.
    @param xt:      Value (probability) of XT events.

    @return times:  List containing the time at which SiPM cells are fired
                    including XT events
    @return idx:    List containing the ID of the hitted cells including
                    XT events.
    """
    if not args.noxt:
        neighbour = [1, -1, CELLSIDE, -CELLSIDE, 1+CELLSIDE, 1-CELLSIDE, -1+CELLSIDE, -1-CELLSIDE]

        for i in range(len(idx)):
            nxt = frandom.randpoiss(xt, 1)[0]
            for j in range(nxt):
                choice = frandom.randint(7, 1)[0]
                idx.append(idx[i] + neighbour[choice])
                times.append(times[i])

    return(times, idx)


def addAP(times, h, ap):
    """!@brief Generation of afterpulses."""
    """!
    Function that generates afterpulses(AP) events and adds them to the list of
    events to simulate. Each avalanche generates a poissonian number of AP
    events. APs are delayed from teyr main signal following a fast and slow
    exponential distributions. Theyr relative signal height is calculated
    consdering the cell recovery time as an RC circuit:
    @f[h=1-e^{-\frac{\Delta_t}{\tau}}@f]

    @param times:   List containing the time at which SiPM cells are fired,
                    including DCR and XT events.
    @param h:       List containing the relative signal height of each
                    fired cell.
    @param ap:      Value (probability) of AP events.

    @return times:  List containing the time at which SiPM cells are fired,
                    including AP events.
    @return h:      List containing the relative signal height of each hitted cell.
    """

    for t in times:
        nap = frandom.randpoiss(ap, 1)[0]
        for j in range(nap):
            delay = frandom.randexp(TAUAPFAST, 1)[0] + frandom.randexp(TAUAPSLOW, 1)[0]
            height = 1 - exp(-delay / CELLRECOVERY)
            times.append(t + delay)
            h.append(height)
    return(times, h)


def SiPMEventAction(times, idx):
    """!@brief Calculates relative signal height for each photoelectron."""
    """!
    Since a SiPM cell may be hitted multiple times, recovery time has to be considered.
    If all cell IDs are different all signal heights are left unchanged with a
    value of 1, else if an ID appears multiple times then the first (in time) hit
    gives a signal height of 1 and the following will have a lower relative height.
    The relative signal height, after a @f$\Delta_t@f$ time from the previous hit,
    is calculated supposing that each SiPM cell recovers as an RC circuit:
    @f[h(\Delta_t)=1-e^{-\frac{t}{\tau}}@f]

    @param times:   List containing the time of all events that have generated
                    an avalanche in the SiPM cells.
    @param idx:     List containing the corresponding cell IDs of the events.
    @return h:      List containing the relative signal height of each avalanche.
    """
    h = [1] * len(times)
    times = np.array(times, dtype='float32', copy=False)
    if not len(idx) == len(set(idx)):
        _, uniqindex, uniqcts = np.unique(idx, return_index=True, return_counts=True)
        for i in range(uniqcts.size):
            if uniqcts[i] > 1:
                midx = idx[uniqindex[i]]
                mtimes = times[idx == midx]
                htemp = 1 - exp(-np.diff(mtimes) / CELLRECOVERY)
                h[i + 1: i + 1 + htemp.size] = htemp
                i += uniqcts[i]
    return h


def SiPMSignalAction(times, sigH, snr=SNR, basespread=0):
    """! @brief Generation of full SiPM signal."""
    """!
    Function that generates the full SiPM signal as the sum of the signals
    of each cell starting from gaussian noise.

    @param times:       List containing the time at wich SiPM cells are fired,
                        including DCR, XT and AP events
    @param sigH:        List containing the correspondin relative pulse height
                        of each fired cell.
    @param snr:         Signal to noise ratio converted into the RMS of the
                        gaussian noise.
    @param basespread:  Sigma of the value to add as baseline spread.
    @return signal:     Array containing the complete sigitized SiPM signal.
    """
    times = np.array(times, dtype='float32', copy=False)
    sigH = np.array(sigH, dtype='float32', copy=False)
    signal = frandom.randn(0, snr, SIGPTS)    # Start with gaussian noise
    if times.size:
        sigH = sigH[times < SIGLEN]
        times = np.uint32(times[times < SIGLEN] / SAMPLING)
        gainvars = frandom.randn(1, CCGV, times.size)    # Each signal has a ccgv
        for i in range(times.size):                      # Generate signals and sum them
            signal += PulseCPU(times[i], sigH[i], gainvars[i])
    return signal

# Fast generation of signals
if args.signal is None:
    ##@cond
    x = np.arange(0, SIGPTS)
    ##@endcond

    ## @var numpy.ndarray $signalmodel
    ## Array containing the signal shape of the SiPM intended as the ideal
    ## signal generated by a signle photoelectron.
    signalmodel = signalgenfortran(0, 1, TFALL / SAMPLING, TRISE / SAMPLING, SIGPTS, 1)

    def PulseCPU(t, h, gainvar):
        """! @brief Generation of single cell signals."""
        """!
        Function that generates the signal from a single SiPM cell.
        This is the "fast" version that uses a pre-computed signal
        shape and moves it in the "right" position in time.

        @param t:       Time at which the cell is triggered.
        @param h:       The relative pulse height of the cell signal
        @param gainvar: Value of cell to cell gain variation for this signal.


        @return s: Array containing the generated cell signal
        """

        # Move the model signal and add ccgv and relative h
        sig = rollfortran(signalmodel, t, gainvar * h)
        return sig


# Full generation of signals (for debug)
else:
    # Generation fully on cpu
    if args.device == 'cpu':
        from libs.libCPU import *
    # Cpu for low light and cpu for high light
    if args.device == 'gpu':
        from libs.libGPU import *


def signalAnalysis(signal, intstart, intgate, threshold):
    """!@brief Features extraction from signals."""
    """!
    Function that extracts some simple features from signals considering an
    integration gate. At the moment the features considered are:
        - Peak: signal peak in the integration gate.
        - Integral: sum of samples in the integration gate corrected by the
        sampling time.
        - Time of Arrival: time position of the first sample above the
        threshold.
        - Time over Threshold: Number fo samples above the threshold corrected
        by the sampling time.
        - Time of Peak: time position of the sample in the peak.
    This list might be expanded in future with more complex features.
    @param signal       Array containing the digitized signal.
    @param intstart     Integer value containing the index corresponding to the
                        start of the integration gate.
    @param intgate      Integer value containing the lenght in units of samples of
                        the integration gate.
    @param threshold    Threshold used to calculate the values described above.
                        It is automatically set to 1.5 times the single
                        photoelectron peak height.
    @return   Returns the features extracted.
    """
    sigingate = signal[intstart:intstart + intgate]

    peak = -1
    integral = -1
    toa = -1
    tot = -1
    top = -1

    if sigingate.max() > threshold:
        peak = sigingate.max()
        integral = sigingate.sum() * SAMPLING
        toa = (sigingate > threshold).argmax() * SAMPLING
        tot = np.count_nonzero(sigingate > threshold) * SAMPLING
        top = (sigingate).argmax() * SAMPLING

    return peak, integral, toa, tot, top


# Non simulation related functions
def somestats(output, realpe=None):
    """!@brief Function that displays histograms of generated events."""
    """!
    @param output:   Array containing the signal features.
    @param realpe:  Optional parameter. Array containing the real number
                        of photoelectrons given to the simulation.

    @sa signalAnalysis for a description of the signal features.
    """

    integral = output[:, 0]
    peak = output[:, 1]
    tstart = output[:, 2]
    tovert = output[:, 3]
    tpeak = output[:, 4]

    plt.figure()
    plt.title('Integral')
    plt.hist(integral, 500, color='k')
    plt.xlabel('Integrated charge [A.U.]')

    plt.figure()
    plt.title('Peak')
    plt.hist(peak, 500, color='k')
    plt.xlabel('Peak value [A.U.]')

    plt.figure()
    plt.title('Starting time')
    plt.hist(tstart, np.arange(tstart.min()-2, tstart.max()+2, 2*SAMPLING), color='k')
    plt.yscale('log')
    plt.xlabel('Starting time [ns]')

    plt.figure()
    plt.title('Time over threshold')
    plt.hist(tovert, np.arange(tovert.min(), tovert.max(), 2*SAMPLING), color='k')
    plt.yscale('log')
    plt.xlabel('Time over threshold [ns]')

    plt.figure()
    plt.title('Peaking time')
    plt.hist(tpeak, np.arange(tpeak.min(), tpeak.max(), 2*SAMPLING), color='k')
    plt.yscale('log')
    plt.xlabel('Peaking time [ns]')

    plt.figure()
    x = np.sort(peak)
    y = np.empty_like(x)
    for i in range(x.size):
        y[i] = np.count_nonzero(x > x[i])

    y *= (1e-3 / y.size / (INTGATE * SAMPLING * 1e-9))
    if peak.max() < 5:
        plt.hlines(DCR/1e3, x.min(), x.max(), 'k', label=f'DCR = {DCR*1e-3:.0f} kHz')
        plt.legend()
    plt.plot(x, y, '.r')
    plt.yscale('log')
    plt.ylabel('Rate [kHz]')
    plt.xlabel('Threshold')
    plt.title('Staircase')

    if realpe is not None:
        plt.figure()
        plt.subplot(121)
        plt.scatter(realpe, peak, c='k', s=2)
        plt.xlabel('Real number of photoelectrons')
        plt.ylabel('Peak value measured')
        plt.subplot(122)
        plt.hist2d(realpe, peak, bins=(100, 100))
        plt.xlabel('Real number of photoelectrons')
        plt.ylabel('Peak value measured')

    plt.figure()
    plt.subplot(121)
    plt.scatter(integral, peak, c='k', s=2)
    plt.xlabel('Peak value measured')
    plt.ylabel('Integrated charge')
    plt.subplot(122)
    plt.hist2d(integral, peak, bins=(100, 100))
    plt.xlabel('Peak value measured')
    plt.ylabel('Integrated charge')
    plt.show()

##@cond
drawn = [False] * nJobs
opened = [True] * nJobs
##@endcond


def sigPlot(signal, npe, ndcr, dev):
    """!Plot each signal pulse produced in the simulation."""
    """!
    Function used for debug purposes only. It creates a window for each
    worker and draws the signals as they are generated.

    @param signal:      Array containing the generated SiPM signal.
    @param npe:         Total number of hitted cells.
    @param ndcr:        Total number of DCR events.
    @param dev:         String that describes the device on which the
                        signal is computed. (cpu / cpu-fast / gpu)
    """

    current_core = multiprocessing.current_process().name
    if current_core == 'MainProcess':
        current_core = 0
    else:
        current_core = int(current_core.split('-')[-1]) - 1

    textstring = f"Core: {current_core:d}\n"
    textstring += f"Device: {dev:s}\n"
    textstring += f"Photons:{npe-ndcr:d} Dcr:{ndcr:d}\n"
    if not drawn[current_core]:
        timearray = np.arange(SIGPTS) * SAMPLING
        ax = plt.subplot(111)
        screenx = 1920
        screeny = 1080
        xsize = int(screenx / 6)
        ysize = int(screeny / 3)
        xpos = (current_core % 6) * xsize
        ypos = ((current_core // 6) % 3) * ysize
        plt.get_current_fig_manager().window.setGeometry(xpos, ypos, xsize, ysize)
        ax.hlines(-0.5, 0, INTSTART * SAMPLING, 'r')
        ax.vlines(INTSTART * SAMPLING, -0.5, -1, 'r')
        ax.vlines((INTSTART - PREG) * SAMPLING, -0.5, -1, 'r')
        ax.hlines(-1, INTSTART * SAMPLING, (INTSTART + INTGATE) * SAMPLING, 'r')
        ax.vlines((INTSTART + INTGATE) * SAMPLING, -1, -0.5, 'r')
        ax.hlines(-0.5, (INTSTART + INTGATE) * SAMPLING, SIGLEN, 'r')
        ax.grid(linestyle=':')
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        drawn[current_core] = True
        ax.plot(timearray, signal, '-b', linewidth=0.5)
    else:
        axlist = plt.gcf().axes
        if len(axlist) == 0:
            opened[current_core] = False
        if opened[current_core]:
            ax = axlist[0]

    if opened[current_core]:
        line = ax.lines[-1]

        txt = ax.text(0.5, 0.70, textstring, transform=ax.transAxes, fontsize=10)

        line.set_ydata(signal)
        ax.relim()
        ax.autoscale_view(tight=False, scalex=False, scaley=True)

        plt.pause(float(args.Graphics)/1000)
        txt.remove()


def initializeRandomPool(seed=None):
    """! Function that initializes random seeds for each worker in the multiprocessing Pool."""
    """!
    This function extracts seeds using os.urandom, so from the operating system
    entropy pool, and uses the to set the seeds for the random generator and
    numpy random generator. Seeds need to be different for each worker,
    otherwise the random values are generate equally among workers."""

    # Get current core number (On MacOs psutil is not working)
    core = multiprocessing.current_process().name
    core = int(core.split('-')[-1])

    # Get some random bits from the sistem entropy pool
    if seed is None:
        seed = int.from_bytes(os.urandom(4), "big")
    seed += core - 1
    # Change rng seed for each worker
    random.seed(seed)
    np.random.seed(seed)
    print("Initializing simulation on worker %s with seed %d\r" % (core, seed))


def SaveFile(fname, out, other=None):
    f = uproot.recreate(fname, compression=uproot.LZ4(9))

    f['SiPMData'] = uproot.newtree({
       'Integral': np.float32,
       'Peak': np.float32,
       'ToA': np.float32,
       'ToT': np.float32,
       'ToP': np.float32})

    f['SiPMData']['Integral'].newbasket(out[:, 0])
    f['SiPMData']['Peak'].newbasket(out[:, 1])
    f['SiPMData']['ToA'].newbasket(out[:, 2])
    f['SiPMData']['ToT'].newbasket(out[:, 3])
    f['SiPMData']['ToP'].newbasket(out[:, 4])

    if other.any():
        other = np.array(other)
        f['GeometryData'] = uproot.newtree({
                'EventId': np.int32,
                'FiberType': np.int8,
                'FiberId': np.int64,
                'FiberX': np.float32,
                'FiberY': np.float32,
                'FiberZ': np.float32
                })

        f['GeometryData']['EventId'].newbasket(np.int32(other[:, 0]))
        f['GeometryData']['FiberType'].newbasket(np.int8(other[:, 1]))
        f['GeometryData']['FiberId'].newbasket(np.int64(other[:, 2]))
        f['GeometryData']['FiberX'].newbasket(np.float32(other[:, 3]))
        f['GeometryData']['FiberY'].newbasket(np.float32(other[:, 4]))
        f['GeometryData']['FiberZ'].newbasket(np.float32(other[:, 5]))


def SaveWaves(fname, signals):
    fname = datetime.now().strftime("%H_%M_%S_") + fname

    sipmsettings = [SIZE,
                    CELLSIZE,
                    DCR,
                    XT,
                    AP,
                    SAMPLING,
                    TRISE * SAMPLING,
                    TFALL * SAMPLING,
                    -20 * np.log10(SNR**2),
                    CCGV]

    with h5py.File(fname, 'a') as hf:
        dset1 = hf.create_dataset('Waveforms',
                                  shape=(signals.shape),
                                  dtype='f',
                                  compression='gzip',
                                  chunks=(1, signals.shape[1]),
                                  compression_opts=9)
        dset2 = hf.create_dataset('SiPMSettings',
                                  shape=(len(sipmsettings),),
                                  dtype='f',
                                  compression='gzip',
                                  compression_opts=9)

        dset1[...] = signals
        dset2[...] = sipmsettings
