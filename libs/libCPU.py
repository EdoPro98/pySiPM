"""In this file I define all the functions I will use in the main file of simulation."""
# EDITING THIS FILE MAY SERIOUSLY COMPROMISE SIMULATION BEHAVIOUR
from libs.FortranFunctions import rollfortran, signalgenfortran, sortfortran
from libs.FortranFunctions import frandom
from variables import *


def PulseCPU(t, h, gainvar):
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
    sig = signalgenfortran(t, h, TFALL / SAMPLING, TRISE / SAMPLING, SIGPTS, gainvar)    # Calculate signal
    return sig
