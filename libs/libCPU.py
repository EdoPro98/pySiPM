# In this file I define all the functions I will use in the main file of simulation
from libs.FortranFunctions import rollfortran, signalgenfortran, sortfortran
from libs.FortranFunctions import frandom
from variables import *

###############################################################################
##>>>   EDITING THIS FILE MAY SERIOUSLY COMPROMISE SIMULATION BEHAVIOUR   <<<##
###############################################################################
def PulseCPU(t, h, gainvar):
    """
    PulseCPU(t,h)

    Function that generates the signal from a single SiPM cell.
    This is the "full" version that computes the signal shape on CPU by evaluating the signal shape function.

    Parameters
    ----------
    t : int32
            Time at which the cell is triggered
    h : float32
            The relative pulse height of the cell signal
    gainvar : float32
            Value of cell to cell gain variation for this signal
    nap : int32
            Number of afterpulses in this signal


    Returns
    -------
    s : np.ndarray
            Array containing the generated cell signal
    """
    sig = signalgenfortran(t, h, TFALL / SAMPLING, TRISE / SAMPLING, SIGPTS, gainvar)    # Calculate signal
    return sig
