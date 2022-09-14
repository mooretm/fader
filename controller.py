"""Controller for simulating DEM gain decreases
    in different frequency bands.
"""

###########
# Imports #
###########
# Import custom modules
from models import Audio
from fader_obj import Fader
import tmsignals as ts


####################
# Initialize Audio #
####################
# Read in audio with SNR of 0
speech_obj = Audio('.\\audio_files_in\\CST_Speech_Trunc.wav', -20)
babble_obj = Audio('.\\audio_files_in\\CST_Babble_4.wav', -20)

# Create signal of interest by combining speech and noise
speech = speech_obj.working_audio
babble = babble_obj.working_audio[0:len(speech)]
combo = speech + babble


#################
# Set constants #
#################
SIGNAL = combo
FS = speech_obj.fs
TRANS_DUR = 3
FLOOR = ts.db2mag(-10)
GAIN = 6
DIRECT_PATH = 'y'


"""Run simulation"""
#######################
# Overall Gain Change #
#######################
oag = Fader(
    signal=SIGNAL,
    fs=FS,
    trans_dur=TRANS_DUR,
    floor=FLOOR,
    gain=GAIN,
    direct_path=DIRECT_PATH
    )

print('-' * 70)
print('OAG')
print('-' * 70)
oag.run(sig_decrease=oag.signal, sig_stable=None)
print('\n')


########################
# Low Freq Gain Change #
########################
lfg = Fader(
    signal=SIGNAL,
    fs=FS,
    trans_dur=TRANS_DUR,
    floor=FLOOR,
    gain=GAIN,
    direct_path=DIRECT_PATH
    )

print('-' * 70)
print('LFG')
print('-' * 70)
lfg.run(sig_decrease=lfg.low, sig_stable=lfg.high)
print('\n')


#########################
# High Freq Gain Change #
#########################
hfg = Fader(
    signal=SIGNAL,
    fs=FS,
    trans_dur=TRANS_DUR,
    floor=FLOOR,
    gain=GAIN,
    direct_path=DIRECT_PATH
    )

print('-' * 70)
print('HFG')
print('-' * 70)
hfg.run(sig_decrease=hfg.high, sig_stable=hfg.low)
print('\n')
