"""Fader utility to change gain over time for selected frequency band.

    Author: Travis M. Moore
    Created: 19 Aug, 2022
    Last Edited: 8 Sep, 2022
"""

###########
# Imports #
###########
# Imports data science packages
import numpy as np
from matplotlib import pyplot as plt

# Import sound packages
import sounddevice as sd
from scipy.io import wavfile

# Import custom modules
from models import Audio
import filters
import tmsignals as ts


#########
# BEGIN #
#########
# Initialize values
# Read in audio and set SNR to +5
speech_obj = Audio('.\\audio_files_in\\CST_Speech_Trunc.wav', -15)
babble_obj = Audio('.\\audio_files_in\\CST_Babble_4.wav', -20)

# Define minimum gain (0 - 1)
FLOOR = 0.325

# Define stable portions
# 0 for transition only (no stable portions)
# 1 for stable onset only
# 2 for stable onset AND offset
STABLE = 2

# The range of delta Ts to use
for t in range(1,21):
    # Set duration parameters in seconds
    stable_dur = 5
    delta_t = t
    total_dur = stable_dur + delta_t + stable_dur

    # Convert seconds to samples
    stable_dur_samps = stable_dur * speech_obj.fs
    delta_t_samps = delta_t * speech_obj.fs
    total_dur_samps = (stable_dur_samps * 2) + delta_t_samps

    # Grab audio clips based on total_dur
    speech = speech_obj.working_audio[0:total_dur_samps]
    babble = babble_obj.working_audio[0:total_dur_samps]
    combo = speech + babble

    plt.plot(combo)
    plt.title('Original Signal')
    #plt.show()

    # Filter audio at 1000 Hz
    combo_low = filters.butter_filt(combo, 'low', 1000, 10, speech_obj.fs, 'n')
    combo_high = filters.butter_filt(combo, 'high', 1000, 10, speech_obj.fs, 'n')

    # Plot signals
    def mk_plot(low, high):
        a1 = plt.subplot2grid((2,2),(0,0), rowspan=1, colspan=1)
        a2 = plt.subplot2grid((2,2),(0,1), rowspan=1, colspan=1)
        a3 = plt.subplot2grid((2,2),(1,0), rowspan=1, colspan=2)
        a1.plot(low, c='b')
        a1.set(title='Freqs < 1kHz', xlabel='Samples', ylabel='Amplitude')
        a1.set_ylim((-1,1))
        a2.plot(high, c='r')
        a2.set(title='Freqs > 1kHz', xlabel='Samples', ylabel='Amplitude')
        a2.set_ylim((-1,1))
        #a3.plot(low, c='b')
        #a3.plot(high, c='r')
        a3.plot(low+high)
        a3.set(title='Recombined Signal', xlabel='Samples', ylabel='Amplitude')
        a3.set_ylim((-1, 1))
        plt.suptitle("Signal Filtering around 1 kHz")
        plt.tight_layout()
        plt.show()

    # Show the high/low-pass filtered signal
    #mk_plot(combo_low, combo_high)


    # Splice audio into start-transition-end clips
    sig_start = combo_low[0:stable_dur_samps]
    sig_trans = combo_low[stable_dur_samps:(stable_dur_samps + delta_t_samps)]
    sig_end = combo_low[(stable_dur_samps+delta_t_samps):((stable_dur_samps*2) 
        + delta_t_samps)]

    # Fade out low-frequencies
    #sig_gated = filters.mk_gate2(combo_low, stable_dur, delta_t, speech_obj.fs)
    start_env = np.ones(len(sig_start))
    trans_env = np.linspace(1, FLOOR, len(sig_trans))
    end_env = np.repeat(FLOOR, len(sig_end))
    
    if STABLE == 2:
        combo_low = combo_low
        envelope = np.hstack([start_env, trans_env, end_env])
    elif STABLE == 1:
        combo_low = combo_low[0:(len(sig_start) + len(sig_trans))]
        envelope = np.hstack([start_env, trans_env])
    elif STABLE == 0:
        combo_low = combo_low[0:len(sig_trans)]
        envelope == trans_env
    
    sig_gated = combo_low * envelope


    # Truncate highs to match filtered lows duration
    combo_high = combo_high[0:len(sig_gated)]

    # Recombine filtered signals
    sig = sig_gated + combo_high

    # Plot 
    #mk_plot(low=sig_gated, high=combo_high)


    # Measure RMS of dry/wet
    # Define regions to measure RMS
    edge1 = (stable_dur) * speech_obj.fs
    edge2 = (stable_dur + delta_t) * speech_obj.fs
    #print(f"edge1: {edge1}")
    #print(f"edge2: {edge2}")

    # Begin plotting signal with polygon highlights
    plt.plot(combo_low, c='b')
    plt.plot(sig_gated, c='g')

    # First RMS region
    coord = [[0,-1], [0,1], [edge1,1], [edge1,-1]]
    coord.append(coord[0]) # repeat the first point to create a 'closed loop'
    # Plot polygon
    xs, ys = zip(*coord) # create lists of x and y values
    plt.fill(xs, ys, edgecolor='none', facecolor="lightsalmon", alpha=0.25)

    if STABLE == 2:
        # Second RMS region
        coord = [[edge2,-1], [edge2,1], [len(sig_gated),1], [len(sig_gated),-1]]
        coord.append(coord[0]) # repeat the first point to create a 'closed loop'
        # Plot polygon
        xs, ys = zip(*coord) # create lists of x and y values
        plt.fill(xs, ys, edgecolor='none', facecolor="lightsalmon", alpha=0.25)

    #plt.show()

    if STABLE == 2:
        start_rms = ts.mag2db(ts.rms(sig_gated[0:edge1]))
        end_rms = ts.mag2db(ts.rms(sig_gated[edge2:]))
        #print(f"Dry RMS: {start_rms}")
        #print(f"Wet RMS: {end_rms}")
        print(f"RMS drop: {abs(start_rms) - abs(end_rms)} dB")

    # Present audio
    #sd.play(sig.T, speech_obj.fs)
    #sd.wait(total_dur)

    wavfile.write('.\\audio_files_out\\LFG_' + str(t) + '.wav', speech_obj.fs, sig)
