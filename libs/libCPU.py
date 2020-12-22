"""In this file I define all the functions I will use in the main file of simulation."""
from variables import *
from libs.FortranFunctions import fsignal


def PulseCPU(t, h):
    """! @brief Generation of single cell signals."""
    """!
    Function that generates the signal from a single SiPM cell.
    This is the "full" version that calculates each signal considering
    a double exponential function.

    @param t:       Time at which the cell is triggered.
    @param h:       The relative pulse height of the cell signal
    @param gainvar: Value of cell to cell gain variation for this signal.


    @return s: Array containing the generated cell signal
    """
    # Calculate signal
    sig = fsignal(t, TF, TR, h*PEAKRATIO, SIGPTS)
    return sig
