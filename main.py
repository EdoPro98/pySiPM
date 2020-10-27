"""!File containing the main function used to simulate SiPM events.

Author: Edoardo Proserpio
Email: eproserpio@studenti.uninsubria.it edoardo.proserpio@gmail.com
"""
# Function of simulation
from libs.lib import *
from libs.FortranFunctions import signalanalysisfortran


def SiPM(times, other=None):
    """
    SiPM(times,other)

    Function that calls all the procedures defined in libs to generate a complete SiPM event.

    Parameters
    -----------
    times : np.ndarray
            This array contains the arriving time of each photon on the SiPM sensor surface.
            This array is the input of the simulation.

    other : tuple
            This variable may contain other informations about the event generated.
            It can be the event id, the arriving time inserted in the simulation
            or the real number of photons inserted in the simulation.
            This tuple will be copied as it is in the output.

    Returns
    ----------
    integral : double
            The integral of the signal calculated in the integration gate

    peak : double
            The height of the signal in the integration gate

    tstart : double
            The time of arrival of the signal in ns defined as the first sample over the threshld of 1.5

    other : tuple
            The copy of the `other` variable given in the input

    signal : np.ndarray(float32)
            If the options -W is enabled the complete SiPM signal will be passed in the output.
            Otherwise this output is None
    """

    # Generate DCR events (times)
    if not args.nodcr:
        dcrTime = addDCR(DCR)
        if dcrTime:
            times.extend(dcrTime)
    # Calculate idx of hitted cells
    idx = HitCells(len(times))
    # Add XT events
    times, idx = addXT(times, idx, XT)
    # Calculate signal height of each cell
    sigH = SiPMEventAction(times, idx)
    times, sigH = addAP(times, sigH, AP)

    # Generate digital signals
    signal = SiPMSignalAction(times, sigH, SNR, BASESPREAD)

    # # Select signal in the integration gate
    integral, peak, tstart, tovert, tpeak = signalAnalysis(signal, INTSTART, INTGATE, THRESHOLD)
    if args.Graphics:
        if not args.signal:
            dev = 'cpu-fast'
        elif args.device == 'cpu':
            dev = 'cpu'
        elif args.device == 'gpu':
            if (len(times) < CPUTHRESHOLD) | (len(times) > GPUMAX):
                dev = 'gpu(cpu)'
            else:
                dev = 'gpu'
        sigPlot(signal, len(times), len(dcrTime), dev)
    if not args.wavedump:
        signal = None
    return(integral, peak, tstart, tovert, tpeak, other, signal)
