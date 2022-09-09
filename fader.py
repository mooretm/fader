"""Fader utility to change gain over time for selected frequency band.

    Author: Travis M. Moore
    Created: 19 Aug, 2022
    Last Edited: 9 Sep, 2022
"""

###########
# Imports #
###########
# Imports data science packages
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

# Import sound packages
import sounddevice as sd
from scipy.io import wavfile

# Import custom modules
from models import Audio
import filters
import tmsignals as ts


#####################
# Initialize values #
#####################
# Read in audio and set SNR to +5
speech_obj = Audio('.\\audio_files_in\\CST_Speech_Trunc.wav', -15)
babble_obj = Audio('.\\audio_files_in\\CST_Babble_4.wav', -20)

# Sampling rate
FS = speech_obj.fs

# Define minimum gain (0 - 1)
FLOOR = 0.325

# Define stable portions
# 0 for transition only (no stable portions)
# 1 for stable onset only
# 2 for stable onset AND offset
STABLE = 2

# The delay for the "amplified" signal
# in relation to the direct pathway signal
# in samples
DELAY = int(np.ceil((5 / 1000) * FS))
print(f"DELAY: {DELAY} samples")


#########
# BEGIN #
#########
# The range of delta Ts to use
for t in range(1,21):
    # Set duration parameters in seconds
    stable_dur = 5
    delta_t = t
    total_dur = stable_dur + delta_t + stable_dur

    # Convert seconds to samples
    stable_dur_samps = stable_dur * FS
    delta_t_samps = delta_t * FS
    total_dur_samps = (stable_dur_samps * 2) + delta_t_samps

    # Grab audio clips based on total_dur
    speech = speech_obj.working_audio[0:total_dur_samps]
    babble = babble_obj.working_audio[0:total_dur_samps]
    combo = speech + babble

    #plt.plot(combo)
    #plt.title('Original Signal')
    #plt.show()

    # Filter audio at 1000 Hz
    combo_low = filters.butter_filt(combo, 'low', 1000, 10, FS, 'n')
    combo_high = filters.butter_filt(combo, 'high', 1000, 10, FS, 'n')
    combo_direct = filters.butter_filt(combo, 'low', 750, 10, FS, 'n')

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
    
    # Create gated signal
    sig_gated = combo_low * envelope

    # Truncate highs to match filtered lows duration
    combo_high = combo_high[0:len(sig_gated)]
    # Truncate direct path sound to match filtered lows duration
    combo_direct = combo_direct[0:len(sig_gated)]

    # Recombine filtered signals
    sig = sig_gated + combo_high

    # Plot 
    #mk_plot(low=sig_gated, high=combo_high)

    # Measure RMS of dry/wet
    # Define regions to measure RMS
    edge1 = (stable_dur) * FS
    edge2 = (stable_dur + delta_t) * FS
    #print(f"edge1: {edge1}")
    #print(f"edge2: {edge2}")


    # Begin plotting signal with polygon highlights
    plt.plot(combo_low, c='b')
    plt.plot(sig_gated, c='g')

    # First RMS region
    if (STABLE == 1) or (STABLE == 2):
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


    ####################
    # Add direct sound #
    ####################
    # Apply delay to sig
    sig = sig[:-DELAY]
    # Truncate direct path to match sig length
    combo_direct = combo_direct[DELAY:]

    # Combine DEM and direct path signals
    final_sig = sig + combo_direct
    final_sig = ts.doNormalize(final_sig, FS)

    # Calculate RMS drop using stable periods
    if STABLE == 2:
        start_rms = ts.mag2db(ts.rms(sig_gated[0:edge1]))
        end_rms = ts.mag2db(ts.rms(sig_gated[edge2:]))
        #print(f"Dry RMS: {start_rms}")
        #print(f"Wet RMS: {end_rms}")
        print(f"DEM RMS drop: {np.round(abs(start_rms) - abs(end_rms), 2)} dB")

        start_rms = ts.mag2db(ts.rms(final_sig[0:edge1]))
        end_rms = ts.mag2db(ts.rms(final_sig[edge2:]))
        #print(f"Dry RMS: {start_rms}")
        #print(f"Wet RMS: {end_rms}")
        print(f"DEM + Direct Path RMS drop: {np.round(abs(start_rms) - abs(end_rms), 2)} dB")


    # Function to plot the delay between DEM and direct path signals
    def mk_plot2():
        plt.subplot(3,1,1)
        plt.plot(sig[0:1000])
        plt.title('Amplified signal')

        plt.subplot(3,1,2)
        plt.plot(combo_direct[0:1000])
        plt.title('Direct pathway signal')

        plt.subplot(3,1,3)
        plt.plot(final_sig[0:1000])
        plt.title('Combined signal')

        plt.show()


    # Function to plot original, DEM and direct path signals
    def mk_plot3():
        plt.plot(combo[0:len(final_sig)], label="Original")
        plt.plot(final_sig, label="DEM + Direct Path")
        plt.plot(sig, label="DEM")
        plt.title('Time Waveforms')
        plt.legend()
        plt.show()


    #mk_plot2()
    #mk_plot3()


    #################
    # Present audio #
    #################
    #sd.play(final_sig.T, FS)
    #sd.wait(total_dur)

    # Write .wav files
    #wavfile.write('.\\audio_files_out\\LFG_' + str(t) + '.wav', speech_obj.fs, sig)
