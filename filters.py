"""Filter-related functions for fader
"""

# Import data science packages
import numpy as np
from scipy import signal


####################
# Filter functions #
####################
def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = signal.butter(order, normal_cutoff, btype='low', analog=False)
    return b, a


def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    #y = signal.lfilter(b, a, data)
    y = signal.filtfilt(b, a, data)
    return y


def butter_filt(sig, type, cutoff, order, fs=48000):
    if type == 'low' or type == 'high':
        nyq = 0.5 * fs
        norm_cutoff = cutoff / nyq
        b, a = signal.butter(order, norm_cutoff, btype=type, analog=False)
        y = signal.filtfilt(b, a, sig)
        return y


def doGate(sig, rampdur=0.02, fs=48000):
    """Apply falling ramp to signal SIG, of 
        duration RAMPDUR.         
    """
    gate =  np.cos(np.linspace(np.pi, 2*np.pi, int(fs*rampdur)))
    # Adjust envelope modulator to be within +/-1
    gate = gate + 1 # translate modulator values to the 0/+2 range
    gate = gate/2 # compress values within 0/+1 range
    gate = np.flip(gate)

    gated = np.array(gate) * np.array(sig)

    return gated
