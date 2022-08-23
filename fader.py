"""Fader utility to change gain over time for selected frequency band.

    Author: Travis M. Moore
    Created: 19 Aug, 2022
    Last Edited: 22 Aug, 2022
"""

###########
# Imports #
###########
# Imports data science packages
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
speech_obj = Audio('CST_Speech.wav', -15)
babble_obj = Audio('CST_Babble_4.wav', -20)

for t in range(1, 20):
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

    # Filter audio
    combo_low = filters.butter_filt(combo, 'low', 1000, 5, speech_obj.fs)
    combo_high = filters.butter_filt(combo, 'high', 1000, 5, speech_obj.fs)

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
        a3.plot(low, c='b')
        a3.plot(high, c='r')
        a3.set(title='Recombined Signal', xlabel='Samples', ylabel='Amplitude')
        a3.set_ylim((-1, 1))
        plt.suptitle("Signal Filtering around 1 kHz")
        plt.tight_layout()
        plt.show()

    #mk_plot(combo_low, combo_high)

    # Splice audio into start-transition-end clips
    sig_start = combo_low[0:stable_dur_samps]
    sig_trans = combo_low[stable_dur_samps:(stable_dur_samps + delta_t_samps)]
    sig_end = combo_low[(stable_dur_samps+delta_t_samps):((stable_dur_samps*2) 
        + delta_t_samps)]

    # Fade out low-frequencies
    #sig_gated = filters.mk_gate(combo_low, total_dur, speech_obj.fs)
    sig_gated = filters.mk_gate2(combo_low, stable_dur, delta_t, speech_obj.fs)

    # Truncate highs to match filtered lows duration
    combo_high = combo_high[0:len(sig_gated)]

    mk_plot(sig_gated, combo_high)

    # Recombine filtered signals
    sig = sig_gated + combo_high


    start_rms = ts.mag2db(ts.rms(sig_gated[0:100000]))
    end_rms = ts.mag2db(ts.rms(sig_gated[550000:]))
    print(f"Dry RMS: {start_rms}")
    print(f"Wet RMS: {end_rms}")


    # Present audio
    #sd.play(sig, speech_obj.fs)
    #sd.wait(total_dur)

    wavfile.write('LFG_' + str(t) + '.wav', speech_obj.fs, sig)
