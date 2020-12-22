"""In this file I define all the functions I will use in the main file of simulation."""
from libs.libCPU import PulseCPU
from libs.FortranFunctions import fsignal
from variables import *

# EDITING THIS FILE MAY SERIOUSLY COMPROMISE SIMULATION BEHAVIOUR

signalShapeGPU = cp.ElementwiseKernel(
    # CUDA kernel that generates the signal.
    # Signals are generated in a matrix, each row is a signal, summed up column-wise
    'int32 x, float32 TFALL, float32 TRISE, float32 h',
    'float32 z',
    'z = h * (__expf(-x / TFALL) - __expf(-x / TRISE));',
    'signalShape')


def PulseGPU(t, h):
    """!@brief Function that generates the signal from all SiPM cells at once."""
    """!
    This is the "full" version that computes the signal shape on GPU. The
    signals are generated at once using a CUDA kernel that generates the
    signals of each cell in a matrix (one cell per row) and then sums everything
    column-wise obtaining the complete SiPM signal.

    @param t:   Array containing times at which each the cell is triggered
    @param h:   Array containing the relative pulse height of each cell signal

    @return s:  Array containing the generated SiPM signal
    """
    # Number of signals to generate
    n = t.size
    # Generate matrix containing times of each fired cell
    vect = (cp.arange(SIGPTS, dtype=cp.int32) + cp.zeros((n, 1), dtype=cp.int32) - t[:, None])
    # Zero before the signal
    vect[vect < 0] = 0
    # Transpose array of height values
    h = h[:, None] * PEAKRATIO
    # Call kernel to generate singal
    sig = cp.sum(signalShapeGPU(vect, TFALL / SAMPLING, TRISE / SAMPLING, h), axis=0)
    # If there are afterpulses generate theyr signals
    return sig


# Function that passes signals times and height to main function for generating signals
def SiPMSignalAction(times, sigH, snr, basespread):
    """! @brief Generation of full SiPM signal."""
    """!
    Function that generates the full SiPM signal as the sum of the signals
    of each cell starting from gaussian noise.
    If there are fewer hitted cells the signal is generated on CPU to avoid
    the overhead caused by moving data to GPU memory. Also if there are too
    many hitted cells the signal is generated on CPU to avoid an out-of-memory
    condition since usually VRAM is smaller than system memory.
    Those conditions are controlled by @ref variables.CPUTHRESHOLD and
    @ref variables.GPUMAX and are set by default at 100 and 2000.

    @param times:       List containing the time at wich SiPM cells are fired,
                        including DCR, XT and AP events
    @param sigH:        List containing the correspondin relative pulse height
                        of each fired cell.
    @param snr:         Signal to noise ratio converted into the RMS of the
                        gaussian noise.
    @param basespread:  Sigma of the value to add as baseline spread.

    @return signal:     Array containing the complete sigitized SiPM signal.
    """
    times = np.array(times, dtype=np.float32, copy=False)
    sigH = np.array(sigH, dtype=np.float32, copy=False)
    sigH = sigH[times < SIGLEN]
    times = (times[times < SIGLEN] / SAMPLING).astype(np.uint32)
    signal = normal(0, SNR, SIGPTS)
    sigH = sigH * normal(1, CCGV, times.size)   # Each signal has a ccgv
    if (times.size < CPUTHRESHOLD) or (times.size > GPUMAX):
        for i in range(times.size):
            signal += PulseCPU(times[i], sigH[i])
        return signal
    else:
        signal = cp.asarray(signal)
        signal += PulseGPU(cp.asarray(times), cp.asarray(sigH))
        return cp.asnumpy(signal)
