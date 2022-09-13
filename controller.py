"""Controller for simulating DEM gain decreases
    in different frequency bands.
"""

###########
# Imports #
###########
# Import custom modules
from models import Audio
from fader_obj import Fader, HFG_Fader, LFG_Fader, OAG_Fader


#########
# BEGIN #
#########
# Read in audio with SNR of 0
speech_obj = Audio('.\\audio_files_in\\CST_Speech_Trunc.wav', -20)
babble_obj = Audio('.\\audio_files_in\\CST_Babble_4.wav', -20)

# Create signal of interest by combining speech and noise
speech = speech_obj.working_audio
babble = babble_obj.working_audio[0:len(speech)]
combo = speech + babble


# OAG
f = OAG_Fader(
    signal=combo, 
    fs=speech_obj.fs,
    trans_dur=1,
    floor=0.25,
    direct_path='y'
)

f.run_sim(f.signal)
#f.plot_segments(sig_ungated=f.signal)
#f.plot_overlay()
#f.play_audio()
f.write_audio()



# LFG
# f = LFG_Fader(
#     signal=combo, 
#     fs=speech_obj.fs,
#     trans_dur=1,
#     floor=0.25,
#     direct_path='y'
# )

# f.run_sim(sig_decrease=f.low, sig_stable=f.high)
# f.plot_high_low(lowpass=f.low, highpass=f.high)
# f.plot_segments(sig_ungated=f.low)
# f.plot_overlay()
# f.play_audio()



# HFG
# f = HFG_Fader(
#     signal=combo, 
#     fs=speech_obj.fs,
#     trans_dur=1,
#     floor=0,
#     direct_path='y'
# )

# f.run_sim(sig_decrease=f.high, sig_stable=f.low)
# f.plot_high_low(lowpass=f.low, highpass=f.high)
# f.plot_segments(sig_ungated=f.high)
# f.plot_overlay()
# f.play_audio()



"""
# Create fader object
fader = Fader(
    signal=combo, 
    fs=speech_obj.fs,
    trans_dur=4,
    floor=0.25,
    direct_path='y'
    )

# Process the input audio
fader.run_sim(sig_decrease=fader.high, sig_stable=fader.low)

# Plots
#fader.mk_plot(lowpass=fader.low, highpass=fader.high)
#fader.mk_plot(lowpass=fader.sig_stable, highpass=fader.sig_gated)
fader.plot_segments()
#fader.plot_delay()
fader.plot_overlay()

# Present audio
#fader.play_audio()

# Write audio
#fader.write_audio()
"""

