#from scipy import interpolate
import numpy as np
from scipy.fft import irfft, rfft, rfftfreq


def addSynth(F0, harm, amp, phi, dur, fs = 48000):
    """ 
        Create a complex signal via additive synthesis. Returns
        the signal AND the time base. 
        
        F0: fundamental frequency in Hz
        HARM: list of F0 harmonics for synthesis
        AMP: amplitude in dB
        PHI: phase in degrees
        DUR: duration in seconds
        FS: sampling rate in samples/second

        EXAMPLE: Create a sawtooth wave
            harms = np.arange(2,60,2)
            amps = [1/x for x in harms]
            phis = np.zeros(len(harms),dtype=int)
            [t, sig] = addSynth(100,harms,amps,phis,0.05,48000)
        
        Written by: Travis M. Moore
        Last edited: Jan. 13, 2022
    """
    phi = deg2rad(phi) # phase to radians
    t = np.arange(0,dur,1/fs) # time base
    sig = np.zeros(len(t)) # Make empty array to save memory
    # Additive synthesis
    for ii in range(len(harm)):
        harmonics = amp[ii] * np.sin(2*np.pi* (F0 * harm[ii]) * t + phi[ii])
        sig += harmonics
    return [t, sig]


def db2mag(db):
    """ 
        Convert decibels to magnitude. Takes a single
        value or a list of values.
    """
    # Must use this form to handle negative db values!
    try:
        mag = [10**(x/20) for x in db]
        return mag
    except:
        mag = 10**(db/20)
        return mag


def deg2rad(deg):
    """ 
        Convert degrees to radians. Takes a single
        value or a list of values.
    """
    try:
        rads = [np.radians(x) for x in deg]
        return rads
    except:
        rads = np.radians(deg)
        return rads


def doFFT(sig,fs,N=2048):
    """
        Wrapper function for scipy rfft. Calculate 
        a single-sided fast Fourier transform for a 
        given signal SIG and sampling rate FS. For 
        a double-sided FFT (i.e., with negative half 
        from imaginary numbers), use scipy fft function.

            SIG: A 1- or 2-channel signal. Note 
                2-channel signals are first combined.
            FS: sampling rate in Hz

        Based on tutorial: https://realpython.com/python-scipy-fft/
        Adapted by: Travis M. Moore
        Last edited: Feb. 2, 2022
    """
    # Combine channels for 2-channel signals
    if len(sig.shape) == 2:
        combo = sig[0] + sig[1]
        sig = combo / combo.max() # normalize

    #N = len(sig)
    yf = rfft(sig, N)
    xf = rfftfreq(N, 1/fs)
    #plt.plot(xf,yf)
    #plt.show()
    return xf, np.abs(yf)


def doGate(sig,rampdur=0.02,fs=48000):
    """
        Apply rising and falling ramps to signal SIG, of 
        duration RAMPDUR. Takes a 1-channel or 2-channel 
        signal. 

            SIG: a 1-channel or 2-channel signal
            RAMPDUR: duration of one side of the gate in 
                seconds
            FS: sampling rate in samples/second

            Example: 
            [t, tone] = mkTone(100,0.4,0,48000)
            gated = doGate(tone,0.01,48000)

        Original code: Anonymous
        Adapted by: Travis M. Moore
        Last edited: Jan. 13, 2022          
    """
    gate =  np.cos(np.linspace(np.pi, 2*np.pi, int(fs*rampdur)))
    # Adjust envelope modulator to be within +/-1
    gate = gate + 1 # translate modulator values to the 0/+2 range
    gate = gate/2 # compress values within 0/+1 range
    # Create offset gate by flipping the array
    offsetgate = np.flip(gate)
    # Check number of channels in signal
    if len(sig.shape) == 1:
        # Create "sustain" portion of envelope
        sustain = np.ones(len(sig)-(2*len(gate)))
        envelope = np.concatenate([gate, sustain, offsetgate])
        gated = envelope * sig
    elif len(sig.shape) == 2:
        # Create "sustain" portion of envelope
        sustain = np.ones(len(sig[0])-(2*len(gate)))
        envelope = np.concatenate([gate, sustain, offsetgate])
        gatedLeft = envelope * sig[0]
        gatedRight = envelope * sig[1]
        gated = np.array([gatedLeft, gatedRight])
    return gated


def doLocatePeak():
    pass







def doLoop(sig,numreps,sildur,fs=48000):
    """ Make a train consisting of NUMREPS of repetitions 
        of a given signal SIG separated by silence (i.e.,
        the ISI) of duration SILDUR. Takes a 1- or 
        2-channel signal. 

            SIG: a 1- or 2-channel signal
            NUMREPS: number of time to repeat SIG
            SILDUR: ISI in seconds

        Written by Travis M. Moore
        Last edited: 1/18/2022    
     """
    if len(sig.shape) == 1:
        # Create silence 
        samps = sildur * fs
        shh = np.zeros(int(samps))
        # Create repetitions
        justone = np.hstack([sig, shh])
        train = np.hstack((justone,) * numreps)
        return train
    elif len(sig.shape) == 2:
        # Create silence 
        samps = sildur * fs
        shh = np.zeros(int(samps))
        shh = np.array([shh,shh]) # Make two channels
        # Create repetitions
        justone = np.hstack([sig, shh])
        train = np.hstack((justone,) * numreps)
        return train


def doNormalize(sig,fs=48000):
    """
        Normalize an array between +1 and -1. 
        Useful for audio signals. Implemented 
        in a way for fast processing because looping
        through values of lengthy signals takes 
        too long for on-the-fly processing.
        Based on: myTarget = [2*(x-np.min(myTarget)) / (np.max(myTarget)-(np.min(myTarget)))-1 for x in myTarget]

            SIG: a 1-channel array
            FS: the sampling rate

        Written by: Travis M. Moore
        Created: May 23, 2022
        Last Edited: May 23, 2022
    """
    # DEPRECATED
    # This method, best case scenario, requires additional 
    # zeroing, and worst case scenario alters (stretches) 
    # the signal to fit between -1 and 1. It's fine to 
    # "normalize" to 1 or -1 and be done. 
    # sig = np.float64(sig) # must specify dtype to avoid overflow
    # #sig = sig-np.min(sig) # have to remove DC offset first
    # denom = np.max(sig) - np.min(sig)
    # sig = sig/denom
    # sig = sig * 2
    # sig = sig -1
    # return sig

    sig = sig - np.mean(sig) # remove DC offset
    sig = sig / np.max(abs(sig)) # normalize
    return sig


    
def mag2db(mag):
    """ 
        Convert magnitude to decibels. Takes a single
        value or a list of values.
    """
    try:
        db = [20 * np.log10(x) for x in mag]
        return db
    except:
        db = 20 * np.log10(mag)
        return db


def mkBinauralNoise(freqs,dur,itd,ild,fs):
#def mkBinauralNoise(freqs,dur,rampdur,itd,ild,fs):
    """ 
        Create FROZEN noise and apply an ITD and/or ILD. Frozen 
        noise will produce the same noise every time the function 
        is called. Update the seed to create new noise. 

            FREQS: list of frequencies to include in noise
            DUR: duration in seconds
            RAMPDUR: rise/fall gate duration in seconds. If 0, 
                then no gating is applied.
            ITD: interaural time difference in microseconds.
                Negative numbers lead to the left. 
            ILD: interaural level difference in dB
                Negative numbers favor the left. 
            FS: sampling rate in samples per second

            EXAMPLE: 
                freqs = np.arange(500,2001)
                sig = ts.mkBinauralNoise(freqs,0.1,0,700,-1,48000)

        Written by: Travis M. Moore
        Last edited: Jan. 13, 2022
    """
    freqs = [x/fs for x in freqs] # Change array into samples
    dur = round(dur * fs) # stim duration in seconds to samples
    itd = fs * itd / 1000000 # unrounded samples

    # Create time bases
    fulldur = np.ceil(dur + np.abs(itd))
    t = np.arange(0,fulldur,dtype=int)
    tlead = t - np.abs(itd)/2 # time vector for leading signal
    tlag = t + np.abs(itd)/2 # time vector for lagging signal

    # Set numpy random seed
    r = np.random.RandomState(12) # only affects local seed; not global
    phi = -2*np.pi + ((2*np.pi)+(2*np.pi))*r.rand(len(freqs)) # using r seed

    # Additive synthesis
    if itd > 0:
        sigLeft = np.zeros(len(tlead))
        for ii in range(len(freqs)):
            harmonics = np.sin(2*np.pi* (1 * freqs[ii]) * tlead + phi[ii])
            sigLeft += harmonics   
        sigRight = np.zeros(len(tlag))
        for ii in range(len(freqs)):
            harmonics = np.sin(2*np.pi* (1 * freqs[ii]) * tlag + phi[ii])
            sigRight += harmonics
    elif itd < 0:
        sigRight = np.zeros(len(tlead))
        for ii in range(len(freqs)):
            harmonics = np.sin(2*np.pi* (1 * freqs[ii]) * tlead + phi[ii])
            sigRight += harmonics   
        sigLeft = np.zeros(len(tlag))
        for ii in range(len(freqs)):
            harmonics = np.sin(2*np.pi* (1 * freqs[ii]) * tlag + phi[ii])
            sigLeft += harmonics
    elif itd == 0:
        sigRight = np.zeros(len(t))
        for ii in range(len(freqs)):
            harmonics = np.sin(2*np.pi* (1 * freqs[ii]) * t + phi[ii])
            sigRight += harmonics   
        sigLeft = np.zeros(len(t))
        for ii in range(len(freqs)):
            harmonics = np.sin(2*np.pi* (1 * freqs[ii]) * t + phi[ii])
            sigLeft += harmonics

    # Apply ILD
    if ild > 0:
        sigLeft = sigLeft / db2mag(np.abs(ild/2))
        sigRight = sigRight * db2mag(np.abs(ild/2))
    elif ild < 0:
        sigLeft = sigLeft * db2mag(np.abs(ild/2))
        sigRight = sigRight / db2mag(np.abs(ild/2))

    # Assign noises to left/right channels
    sigBoth = np.array([sigLeft, sigRight])

    """
    # Apply gates
    if rampdur == 0:
        pass
    else:
        # Ensure rampdur is in SECONDS, not samples
        sigBoth = doGate(sigBoth,rampdur)
    """

    return sigBoth


def mkGaborClick(cf, dur, itd, ild, fs):
    """ 
        MKGABORCLICKS Generate two Gabor clicks with a 
        specified ITD/ILD

        CF = carrier frequency in Herz
        DUR = clicklength in seconds
        ITD = delay between clicks in microseconds
        ILD = level difference between L/R channels
        FS = sampling rate in samples/second

        Example: gclick = mkGaborClick(4000,0.002,300,2,48000);

       Author: Chris Stecker
       Adapted by: Travis Moore
       Date of MATLAB adaptation: Apr. 19, 2017
       Date of Python adaptation: Jan. 11, 2022
       Last Edited: Jan. 11, 2022
    """
    # Convert user-friendly argument values to samples:
    cf = cf / fs
    dur = dur * fs
    itd = fs * itd / 1000000
    sd = dur / 6.66; # samples per sd. 6.66 sd total dur (full 16 bits)

    # Is this conversion necessary?
    #ILD = 20*np.log10(db2mag(abs(ild)))

    # Time base
    fulldur = dur + np.abs(itd)
    t = np.arange(1,fulldur+1) - (fulldur/2) # have to add 1 due to arange not include final value
    tlead = t - itd/2 # assuming t and itd in samps
    tlag = t + itd/2 # assuming t and itd in samps

    # Create clicks
    clickleft = np.cos(2*np.pi*cf*tlead) * np.exp(np.square((tlead/sd))*-1)
    clickright = np.cos(2*np.pi*cf*tlag) * np.exp(np.square((tlag/sd))*-1)

    """ See newer version immediately below 
    # Apply ILD
    if ild > 0:
        clickleft = clickleft / ILD
    if ild < 0:
        clickright = clickright / ILD
    """
    # Apply ILD
    if ild > 0:
        clickleft = clickleft / db2mag(np.abs(ild/2))
        clickright = clickright * db2mag(np.abs(ild/2))
    elif ild < 0:
        clickleft = clickleft * db2mag(np.abs(ild/2))
        clickright = clickright / db2mag(np.abs(ild/2))

    clickBoth = np.array([clickleft, clickright])
    return clickBoth


def mkIPD(freq,dur,ipd,ild,fs=48000):
    """
        Create a binaural pure tone at frequency FREQ 
        with an interaural phase delay (IPD) and/or interaural 
        level difference (ILD). 

            FREQ: frequency in Hz
            DUR: duration in seconds
            GATEDUR: duration of a single ramp in seconds
            IPD: phase difference in DEGREES
            ILD: level difference in dB

            EXAMPLE: sig = mkIPD(500,0.1,90,-2,48000)

        Written by: Travis M. Moore
        Last edited: Jan. 19, 2022
    """
    # Apply +/- half of phase to each channel
    phiHalf = np.abs(ipd)/2
    if ipd <= 0:
        [t, lchan] = mkTone(freq,dur,np.abs(phiHalf),fs) # wants deg
        [t, rchan] = mkTone(freq,dur,(phiHalf*-1),fs) # wants deg
    elif ipd > 0:
        [t, lchan] = mkTone(freq,dur,(phiHalf*-1),fs) # wants deg
        [t, rchan] = mkTone(freq,dur,phiHalf,fs) # wants deg

    # Apply ILD
    if ild > 0:
        lchan = lchan / db2mag(np.abs(ild/2))
        rchan = rchan * db2mag(np.abs(ild/2))
    elif ild < 0:
        lchan = lchan * db2mag(np.abs(ild/2))
        rchan = rchan / db2mag(np.abs(ild/2))

    # Assign channels
    sig2chan = np.array([lchan, rchan])
    return sig2chan


def mkITD(freq,dur,itd,ild,rampdur,fs=48000):
    """
        Create a binaural pure tone at frequency FREQ with 
        an interaural time delay (ITD) and/or interaural 
        level difference (ILD). Implements a whole-waveform
        shift in time. 

            FREQ: frequency in Hz
            DUR: duration in seconds
            ITD: time difference in MICROSECONDS
            ILD: level difference in dB
            RAMPDUR: duration of gating in seconds
            FS: sampling rate in Hz

            EXAMPLE: sig = mkITD(500,0.05,800,0,0.02,48000)

        Written by: Travis M. Moore
        Last edited: Feb. 8, 2022
    """
    freq = freq / fs # CPS to cycles per sample
    dur = round(dur * fs) # stim duration in seconds to samples
    itd = fs * itd / 1000000 # unrounded samples
    itd_int = int(np.ceil(np.abs(itd)/2)) # rounded itd samples
    rampdur = round(rampdur * fs)

    # Create time vectors
    fulldur = np.ceil(dur + np.abs(itd))
    t = np.arange(0,fulldur,dtype=int)
    tlead = t - (np.abs(itd)/2) # time vector for leading signal
    tlag = t + (np.abs(itd)/2) # time vector for lagging signal

    #gate = np.cos(np.linspace(np.pi, 2*np.pi, int(fs*rampdur)))
    gate = np.cos(np.linspace(np.pi, 2*np.pi, rampdur))
    # Adjust envelope modulator to be within +/-1
    gate = gate + 1 # translate modulator values to the 0/+2 range
    gate = gate/2 # compress values within 0/+1 range
    # Create offset gate by flipping the array
    offsetgate = np.flip(gate)
    # Create "sustain" portion of envelope
    sustain = np.ones(len(tlead)-(2*len(gate))-itd_int)
    # Create padding for itd offset
    pad = np.zeros(itd_int)

    # Create gated ITD signal
    if itd < 0:
        envelopeLeft = np.concatenate([gate, sustain, offsetgate, pad])
        lchan = np.sin(2*np.pi*freq*tlag) * envelopeLeft
        envelopeRight = np.concatenate([pad, gate, sustain, offsetgate])
        rchan = np.sin(2*np.pi*freq*tlead) * envelopeRight
    elif itd > 0:
        envelopeLeft = np.concatenate([pad,gate, sustain, offsetgate])
        lchan = np.sin(2*np.pi*freq*tlead) * envelopeLeft
        envelopeRight = np.concatenate([gate, sustain, offsetgate, pad])
        rchan = np.sin(2*np.pi*freq*tlag) * envelopeRight
    elif itd == 0:
        envelope = np.concatenate([gate, sustain, offsetgate])
        lchan = np.sin(2*np.pi*freq*t) * envelope
        rchan = np.sin(2*np.pi*freq*t) * envelope
    
    # Apply ILD
    if ild < 0:
        lchan = lchan * db2mag(np.abs(ild/2))
        rchan = rchan / db2mag(np.abs(ild/2))
    elif ild > 0:
        lchan = lchan / db2mag(np.abs(ild/2))
        rchan = rchan * db2mag(np.abs(ild/2))

    sigBoth = np.array([lchan, rchan])
    return sigBoth


def mkNoise(freqs,dur,fs):
    """ Create a brief noise using additive synthesis and random 
        phases. WARNING: Durations longer than 1 second will 
        repeat. Large bandwidths can take a while to process.
        
            FREQS: a list of frequencies to be included in the noise
            DUR: duration in seconds
            FS: sampling rate in Hz

            EXAMPLE: sig = mkNoise(np.arange(250,3010,10),0.5,48000)

        Written by Travis M. Moore
        Last edited: 1/13/2022
    """
    # Create arbitrary amplitudes for addSynth
    amp = np.ones(len(freqs))
    # Create random phases
    phi = -2*np.pi + ((2*np.pi)+(2*np.pi))*np.random.rand(len(freqs))
    phi = rad2deg(phi) # addSynth function wants degrees
    [t, myNoise] = addSynth(1,freqs,amp,phi,dur,fs)
    return myNoise


def mkTone(freq, dur, phi=0, fs=48000):
    """ Create a pure tone. Returns the signal 
        AND the time base. 
    
        FREQ: frequency in Hz
        DUR: duration in SECONDS
        PHI: phase in DEGREES
        FS: sampling rate

        EXAMPLE: [t, sig] = (500,0.2,0,48000)

    Written by: Travis M. Moore
    Last edited: 1/12/2022
    """
    phi = np.deg2rad(phi) # to radians
    t = np.arange(0,dur,1/fs) # time base
    sig = np.sin(2*np.pi*freq*t+phi)
    return [t, sig]


def phase2time(freq,phase):
    """
        Calculate a time difference in MICROSECONDS from 
        a given phase in DEGREES. 

            FREQ: frequency in Hz
            PHASE: the phase in DEGREES

        Written by: Travis M. Moore
        Last edited: Jan. 13, 2022
    """
    itd = phase / (360*freq)
    return itd*1000000


def rad2deg(rad):
    """ Convert radians to degrees. Takes a single
        value or a list of values.
    
    Written by: Travis M. Moore
    Last edited: 1/12/2022
    """
    try:
        degs = [np.degrees(x) for x in rad]
        return degs
    except:
        degs = np.degrees(rad)
        return degs


def rms(sig):
    """ 
        Calculate the root mean square of a signal. 
        
        NOTE: np.square will return invalid, negative 
            results if the number excedes the bit 
            depth. In these cases, convert to int64
            EXAMPLE: sig = np.array(sig,dtype=int)

        Written by: Travis M. Moore
        Last edited: Feb. 3, 2020
    """
    theRMS = np.sqrt(np.mean(np.square(sig)))
    return theRMS


def setRMS(sig,amp,eq='n'):
    """
        Set RMS level of a 1-channel or 2-channel signal.
    
        SIG: a 1-channel or 2-channel signal
        AMP: the desired amplitude to be applied to 
            each channel. Note this will be the RMS 
            per channel, not the total of both channels.
        EQ: takes 'y' or 'n'. Whether or not to equalize 
            the levels in a 2-channel signal. For example, 
            a signal with an ILD would lose the ILD with 
            EQ='y', so the default in 'n'.

        EXAMPLE: 
        Create a 2 channel signal
        [t, tone1] = mkTone(200,0.1,30,48000)
        [t, tone2] = mkTone(100,0.1,0,48000)
        combo = np.array([tone1, tone2])
        adjusted = setRMS(combo,-15)

        Written by: Travis M. Moore
        Created: Jan. 10, 2022
        Last edited: May 17, 2022
    """
    if len(sig.shape) == 1:
        rmsdb = mag2db(rms(sig))
        refdb = amp
        diffdb = np.abs(rmsdb - refdb)
        if rmsdb > refdb:
            sigAdj = sig / db2mag(diffdb)
        elif rmsdb < refdb:
            sigAdj = sig * db2mag(diffdb)
        # Edit 5/17/22
        # Added handling for when rmsdb == refdb
        elif rmsdb == refdb:
            sigAdj = sig
        return sigAdj
        
    elif len(sig.shape) == 2:
        rmsdbLeft = mag2db(rms(sig[0]))
        rmsdbRight = mag2db(rms(sig[1]))

        ILD = np.abs(rmsdbLeft - rmsdbRight) # get lvl diff

        # Determine lvl advantage
        if rmsdbLeft > rmsdbRight:
            lvlAdv = 'left'
            #print("Adv: %s" % lvlAdv)
        elif rmsdbRight > rmsdbLeft:
            lvlAdv = 'right'
            #print("Adv: %s" % lvlAdv)
        elif rmsdbLeft == rmsdbRight:
            lvlAdv = None

        #refdb = amp - 3 # apply half amp to each channel
        refdb = amp
        diffdbLeft = np.abs(rmsdbLeft - refdb)
        diffdbRight = np.abs(rmsdbRight - refdb)

        # Adjust left channel
        if rmsdbLeft > refdb:
            sigAdjLeft = sig[0] / db2mag(diffdbLeft)
        elif rmsdbLeft < refdb:
            sigAdjLeft = sig[0] * db2mag(diffdbLeft)
        # Adjust right channel
        if rmsdbRight > refdb:
            sigAdjRight = sig[1] / db2mag(diffdbRight)
        elif rmsdbRight < refdb:
            sigAdjRight = sig[1] * db2mag(diffdbRight)

        # If there is a lvl difference to maintain across channels
        if eq == 'n':
            if lvlAdv == 'left':
                sigAdjLeft = sigAdjLeft * db2mag(ILD/2)
                sigAdjRight = sigAdjRight / db2mag(ILD/2)
            elif lvlAdv == 'right':
                sigAdjLeft = sigAdjLeft / db2mag(ILD/2)
                sigAdjRight = sigAdjRight * db2mag(ILD/2)

        sigBothAdj = np.array([sigAdjLeft, sigAdjRight])
        return sigBothAdj


def specLvl(sig, upr, lwr):
    """
        Calculate the spectral density level of a given
        noise. 

            SIG: a signal consisting of a bandwidth > 1 Hz
            UPR: the upper bandwidth cutoff frequency
            LWR: the lower bandwidth cutoff frequency

        Written by: Travis M. Moore
        Last edited: 2019
    """
    OAL = mag2db(rms(sig))
    BW = upr - lwr
    SPLnb = OAL - np.log10(BW/1)
    return SPLnb


def time2phase(freq,itd):
    """
        Calculate a phase shift in DEGREES from an ITD
        in MICROSECONDS. 

            FREQ: frequency in Hz
            ITD: the time difference in MICROSECONDS

        Written by: Travis M. Moore
        Last edited: Jan. 13, 2022
    """
    ps = 360 * itd * freq
    return ps/1000000
