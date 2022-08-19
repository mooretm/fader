"""Filter-related functions for fader
"""

# Import data science packages
import numpy as np
from scipy import signal


####################
# Filter functions #
####################
def butter_filt(sig, type, cutoff, order, fs=48000):
    nyq = 0.5 * fs
    norm_cutoff = cutoff / nyq
    b, a = signal.butter(order, norm_cutoff, btype=type, analog=False)
    y = signal.filtfilt(b, a, sig)
    return y


def mk_gate(sig, win_dur, fs):
    """Apply falling ramp to signal SIG, of 
        duration RAMPDUR.         
    """
    gate =  np.cos(np.linspace(np.pi, 2*np.pi, int(fs*win_dur)))
    # Adjust envelope modulator to be within +/-1
    gate = gate + 1 # translate modulator values to the 0/+2 range
    gate = gate/2 # compress values within 0/+1 range
    gate = np.flip(gate)

    gated = np.array(gate) * np.array(sig)

    return gated
