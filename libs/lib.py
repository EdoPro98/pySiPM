# In this file I define all the functions I will use in the main file of simulation
from libs.FortranFunctions import rollfortran, signalgenfortran, sortfortran
from libs.FortranFunctions import frandom
from variables import *


###############################################################################
##>>>   EDITING THIS FILE MAY SERIOUSLY COMPROMISE SIMULATION BEHAVIOUR   <<<##
###############################################################################

### GENERATION OF DCR TIMES ###
def addDCR(rate):
    """
    addDCR(rate)

    Function that generates times of dcr events.

    Parameters
    ------------
    rate : double
            Rate of dcr in kHz

    Returns
    --------
    dcrTime : np.ndarray
            Array containing dcr times
    """
    dcrTime = [0]
    while dcrTime[-1] < SIGLEN:  # Generate from exponential distribution of delays
        delayDCR = frandom.randexp(1/rate, 1) * 1e9
        dcrTime.extend(dcrTime[-1] + delayDCR)

    dcrTime = np.array(dcrTime[1:-1])
    return dcrTime


### FUNCTION TO UPDATE THE SIPM MATRIX ARRAY ###
def HitCells(times):
    """
    HitCells(times)

    Function that generates the cell IDs for the corresponding times.

    Parameters
    ----------
    time : np.ndarray
        Array containing the time at which SiPM cells are fired

    Returns
    -------
    idx : np.ndarray(int16)
            Array containing the ID of the hitted cells
    """
    idx = frandom.randint(NCELL, times.size)
    return(np.uint16(idx))


def addXT(times, idx, xt):
    """
    addXT(times, idx, xt)

    Function that generates the times and cell IDs for the XT events and adds them to the list.

    Parameters
    ----------
    times : np.ndarray
        Array containing the time at which SiPM cells are fired
    idx : np.ndarray
        Array containing the ID of each fired cell
    xt : double
        Value (probability) of XT events

    Returns
    -------
    times : np.ndarray
        Array containing the time at which SiPM cells are fired
    idx : np.ndarray(int16)
            Array containing the ID of the hitted cells
    """
    if not args.noxt:
        niter = 3
        xtgen = list(idx)
        xtgent = list(times)
        neighbour = [1, -1, CELLSIDE, -CELLSIDE, 1+CELLSIDE, 1-CELLSIDE, -1+CELLSIDE, -1-CELLSIDE]

        for i in range(niter):
            temp = []
            tempt = []
            for j in range(len(xtgen)):
                nxt = frandom.randpoiss(xt, 1)[0]
                for k in range(nxt):
                    choose = frandom.randint(7, 1)[0]
                    temp.append(xtgen[j] + neighbour[choose])
                    tempt.append(xtgent[j])
            if len(temp) == 0:
                break
            xtgen = temp
            xtgent = tempt
            idx = hstack((idx, temp))
            times = hstack((times, tempt))
    return(times, idx)


def addAP(times, h, ap):
    """
    addAP(times, idx, ap)

    Function that generates the times and heights for the AP events and adds them to the list.

    Parameters
    ----------
    times : np.ndarray
        Array containing the time at which SiPM cells are fired
    h : np.ndarray
        Array containing the signal height of each fired cell
    ap : double
        Value (probability) of AP events

    Returns
    -------
    times : np.ndarray
        Array containing the time at which SiPM cells are fired
    h : np.ndarray
            Array containing the signal height of the hitted cells
    """
    temp = []
    temph = []
    nap = frandom.randpoiss(ap, times.size)
    for i, t in enumerate(times):
        for j in range(nap[i]):
            delay = frandom.randexp(TAUAPFAST, 1)[0] + frandom.randexp(TAUAPSLOW, 1)[0]
            height = 1 - exp(-delay / CELLRECOVERY)
            temp.append(t + delay)
            temph.append(height)
    if len(temp):
        times = hstack((times, temp))
        h = hstack((h, temph))
    return(times, h)


def SiPMEventAction(times, idx):
    h = times * 0 + 1
    if not idx.size == len(set(idx)):
        _, uniqindex, uniqcts = np.unique(idx, return_index=True, return_counts=True)
        for i in range(uniqcts.size):
            if uniqcts[i] > 1:
                midx = idx[uniqindex[i]]
                mtimes = times[idx == midx]
                htemp = 1 - exp(-np.diff(mtimes) / CELLRECOVERY)
                h[i] = 1
                h[i + 1: i + 1 + htemp.size] = htemp
                i += uniqcts[i]
    return h


### GENERATION OF SIGNALS SHAPES FAST ###
def SiPMSignalAction(times, sigH, SNR, BASESPREAD):
    """
    signalGen(times,sigH,SNR,BASESPREAD)

    Function that passes signal height and times to the main function that
    generates single signals. Also adds noise.

    Parameters
    ----------
    times : np.ndarray(int32)
    Array containing the time at wich SiPM cells are fired, including xt events (sorted)
    sigH : np.ndarray(float32)
    Array containing the pulse height of each fired SiPM cell
    SNR : double
    The signal to noise ratio of the noise to add
    BASESPREAD : double
    Sigma of the value to add as baseline

    Returns
    --------
    signal : np.ndarray
    Array containing the generated SiPM signal
    """
    baseline = 0
    signal = frandom.randn(baseline, SNR, SIGPTS)    # Start with gaussian noise
    if times.size:
        sigH = sigH[times < SIGLEN]
        times = np.uint32(times[times < SIGLEN] / SAMPLING)
        gainvars = frandom.randn(1, CCGV, times.size)    # Each signal has a ccgv
        for i in range(times.size):                      # Generate signals and sum them
            signal += PulseCPU(times[i], sigH[i], gainvars[i])
    return signal


if args.signal is None:  # If generating signals fast (default)
    x = np.arange(0, SIGPTS)
    # Define the model of my signal (calculate it only once)
    signalmodel = signalgenfortran(0, 1, TFALL, TRISE, SIGPTS, 1)


    def PulseCPU(t, h, gainvar):
        """
        PulseCPU(t,h)

        Function that generates the signal from a single SiPM cell.
        This is the "fast" version that uses a pre-computed signal
        shape and translates it in time.

        Parameters
        ----------
        t : int32
                Time at which the cell is triggered
        h : float32
                The relative pulse height of the cell signal
        gainvar : float32
                Value of cell to cell gain variation for this signal


        Returns
        -------
        s : np.ndarray
                Array containing the generated cell signal
        """

        sig = rollfortran(signalmodel, t, gainvar, h)   # Move the model signal
        return sig


### FULL GENERATION OF SIGNAL SHAPES ###
else:								# Recalculate the signal each time
    if args.device == 'cpu':
        from libs.libCPU import *   # Generation fully on cpu
    if args.device == 'gpu':
        from libs.libGPU import *   # Cpu for low light and cpu for high light


def signalAnalysis(signal, intstart, intgate, threshld):
    sigingate = signal[intstart:intstart + intgate]

    peak = sigingate.max()
    if peak > threshld:
        integral = sigingate.sum() * SAMPLING
        toa = (sigingate > threshld).argmax() * SAMPLING
        tot = np.count_nonzero(sigingate > threshld) * SAMPLING
        top = (sigingate).argmax() * SAMPLING
    else:
        peak = -1
        integral = -1
        toa = -1
        tot = -1
        top = -1
    return peak, integral, toa, tot, top



### SOME STATISCTICS AT END OF SCRIPT ###
def somestats(output, realpe=None):
    """
    somestats(output)

    Function that displays histograms of generated events

    Parameters
    ----------
    output : np.ndarray
            Array containing output of the simulation.
            This array contains the integral,
            peak and starting time in the first three columns.

    Note
    -----
    See docs for a detailed description of integral, peak and tstart
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
    plt.xlabel('Starting time [ns]')
    plt.yscale('log')

    plt.figure()
    plt.title('Time over threshold')
    plt.hist(tovert, np.arange(tovert.min(), tovert.max(), 2*SAMPLING), color='k')
    plt.xlabel('Time over threshold [ns]')

    plt.figure()
    plt.title('Peaking time')
    plt.hist(tpeak, np.arange(tpeak.min(), tpeak.max(), 2*SAMPLING), color='k')
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
    plt.scatter(peak, integral, c='k', s=2)
    plt.xlabel('Peak value measured')
    plt.ylabel('Integrated charge')
    plt.subplot(122)
    plt.hist2d(peak, integral, bins=(100, 100))
    plt.xlabel('Peak value measured')
    plt.ylabel('Integrated charge')
    plt.show()


drawn = [False] * nJobs
opened = [True] * nJobs
def sigPlot(signal, sigTimes, dcrTime, dev):
    """
    sigPlot(signal,sigTimes,dcrTime,dev)

    Function that plots each signal pulse produced in the simulation.

    Parameters
    -----------
    signal : np.ndarray
            Array containing the generated SiPM signal
    sigTimes : np.ndarray
            Array containing all photons events times (including xt)
    dcrTime : np.ndarray
            Array containing dcr events times
    dev : str
            String that describes the device on which the signal is computed
    """

    current_core = multiprocessing.current_process().name
    if current_core == 'MainProcess':
        current_core = 0
    else:
        current_core = int(current_core.split('-')[-1]) - 1

    textstring = f"Core: {current_core:d}\n"
    textstring += f"Device: {dev:s}\n"
    textstring += f"Photons:{sigTimes.size-dcrTime.size:d} Dcr:{dcrTime.size:d}\n"
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


def initializeRandomPool():
    """
    initializeRandomPool()

    Function that initializes random seeds for each worker in the multiprocessing Pool
    """

    # Get current core number (On MacOs psutil is not working)
    current_core = multiprocessing.current_process().name
    current_core = int(current_core.split('-')[-1])
    # Get some random bits from the sistem entropy pool
    rngseed = struct.unpack('I', os.urandom(4))[0] + current_core
    random.seed(rngseed)   # Change rng seed for each worker
    np.random.seed(rngseed)
    core = multiprocessing.current_process().name
    print("Initializing simulation on %s with seed %d\r" % (core, rngseed))


def SaveFile(fname, out, other=None):
    f = uproot.recreate(fname, compression=uproot.LZ4(5))

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
    while os.path.exists(fname):
        f = fname.split('.')
        if f[0][-1].isnumeric():
            fname = f[0][:-1]+str(int(f[0][-1])+1)+'.hdf5'
        else:
            fname = f[0]+'0'+'.'+f[1]

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
