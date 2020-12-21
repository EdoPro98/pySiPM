"""!File containing the main function used to simulate SiPM events.

Author: Edoardo Proserpio
Email: eproserpio@studenti.uninsubria.it edoardo.proserpio@gmail.com
"""
# Function of simulation
from libs.lib import *
from libs.FortranFunctions import signalanalysisfortran


def SiPM(times, other=None):
    """! @biref that calls all the procedures defined in libs to generate a complete SiPM event."""
    """! This function is the main function used to simulate a complete SiPM
     event. It calls all the other methods with the correct parameters and user
     settings.

    @param times This list contains the arriving time of each photon on the SiPM sensor surface.
    This list is the main input of the simulation.

    @param other This variable may contain other informations about the event
    generated. It can be the event id, the arriving time inserted in the
    simulation or the real number of photons inserted in the simulation.
    This tuple will be copied as it is in the output.

    @return integral The integral of the signal calculated in the integration
     gate.

    @return peak The height of the signal in the integration gate.

    @return tstart The time of arrival of the signal in ns defined as the first
     sample over the threshld of 1.5

    @return other The copy of the `other` variable given in the input.

    @return signal If the options -W is enabled the complete SiPM signal
    will be passed in the output. Otherwise this output is "None".
    """
    npe = len(times)
    ndcr = 0
    nxt = 0
    nap = 0

    # Generate DCR events (times)
    if not args.nodcr:
        dcrTime, ndcr = addDCR(DCR)
        if ndcr:
            times.extend(dcrTime)

    # Calculate idx of hitted cells
    idx = HitCells(len(times))

    # Add XT events
    if not args.noxt:
        times, idx, nxt = addXT(times, idx, XT)

    # Calculate signal height of each cell
    times, sigH = SiPMEventAction(times, idx)

    # Add AP events
    if not args.noap:
        times, sigH, nap = addAP(times, sigH, AP)

    # Generate digital signals
    signal = SiPMSignalAction(times, sigH, SNR, BASESPREAD)

    # # Select signal in the integration gate
    peak, integral, tstart, tovert, tpeak = signalAnalysis(signal, INTSTART, INTGATE, THRESHOLD)

    # Plots
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
        sigPlot(signal, len(times), ndcr, dev)

    debug = (npe, ndcr, nxt, nap)
    return (integral, peak, tstart, tovert, tpeak), other, signal, debug
