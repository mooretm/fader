"""Fader utility to change gain over time for selected frequency band.

    Author: Travis M. Moore
    Created: 19 Aug, 2022
    Last Edited: 19 Aug, 2022
"""

# Imports data science packages
import numpy as np
from matplotlib import pyplot as plt

# Import sound packages
import sounddevice as sd

# Import custom modules
from models import Audio
from filters import doGate
from filters import butter_lowpass_filter

speech = Audio('CST_Speech.wav', -15)
babble = Audio('CST_Babble_1.wav', -20)


stable_dur = 5
delta_t = 10
total_dur = stable_dur + delta_t


stable_dur_samps = stable_dur * speech.fs
speech_temp = speech.working_audio[0:stable_dur_samps]
babble_temp = babble.working_audio[0:stable_dur_samps]
combo_temp = speech_temp + babble_temp

combo_gated = doGate(combo_temp, stable_dur, fs=speech.fs)

sd.play(combo_gated, speech.fs)
sd.wait(stable_dur)






