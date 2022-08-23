"""Filter-related functions for fader
"""

# Import data science packages
import numpy as np
from scipy import signal

# Import custom modules
from tmsignals import deg2rad, rad2deg, mag2db, db2mag


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


def mk_gate2(sig, steady_dur, trans_dur, fs):
    """Apply falling ramp to signal SIG, of 
        duration RAMPDUR.         
    """
    # Ungated steady
    stable1 = np.ones(steady_dur*fs, dtype=int)

    # Transition
    #gate =  np.cos(np.linspace(np.pi, 2*np.pi, int(fs*trans_dur)))
    gate = np.cos(np.linspace(deg2rad(360), deg2rad(240), int(fs*trans_dur)))


    # Adjust envelope modulator to be within +/-1
    gate = gate + 1 # translate modulator values to the 0/+2 range
    gate = gate/2 # compress values within 0/+1 range
    #gate = np.flip(gate)

    # Gated steady
    stable2 = np.repeat(gate[-1], steady_dur*fs)

    envelope = np.concatenate([stable1, gate, stable2])

    sig = sig[0:len(envelope)]

    gated = np.array(envelope) * np.array(sig)

    return gated
