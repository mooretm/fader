"""Classes for fader
"""

###########
# Imports #
###########
# Import data science packages
import numpy as np
from scipy.io import wavfile

# Import system packages
import os

# Import sound packages
import sounddevice as sd


#################
# BEGIN CLASSES #
#################
class Audio:
    """ An object for use with .wav files. Audio objects 
        can read a given .wav file, handle audio data type 
        conversion, and store information about a .wav 
        file.
    """
    # Dictionary of data types and ranges for conversions
    wav_dict = {
        'float32': (-1.0, 1.0),
        'int32': (-2147483648, 2147483647),
        'int16': (-32768, 32767),
        'uint8': (0, 255)
    }

    def __init__(self, file_path, level):
        # Parse file path
        self.directory = file_path.split(os.sep) # path only
        self.name = str(file_path.split(os.sep)[-1]) # file name only
        self.file_path = file_path
        self.level = level

        # Read audio file
        fs, audio_file = wavfile.read(self.file_path)

        # Assign audio file attributes
        self.fs = fs
        self.original_audio = audio_file
        self.dur = len(self.original_audio) / self.fs
        self.t = np.arange(0,self.dur, 1/self.fs)
        self.data_type = np.dtype(audio_file[0])
        print(f"Incoming audio data type: {self.data_type}")

        # Immediately convert to float64 for processing
        self.convert_to_float()


    def convert_to_float(self):
        """ Convert original audio data type to float64 
            for processing
        """
        if self.data_type == 'float64':
            pass
        else:
            # 1. Convert to float64
            sig = self.original_audio.astype(np.float64)
            # 2. Divide by original dtype max val
            sig = sig / self.wav_dict[str(self.data_type)][1]
            self.working_audio = sig


    def play(self):
        """ Present working audio """
        print(f"Presenting audio data type: {np.dtype(self.working_audio[0])}")
        # plt.subplot(1,3,1)
        # plt.plot(self.original_audio)
        # plt.subplot(1,3,2)
        # plt.plot(self.working_audio)

        sig = self.setRMS()
        self.working_audio = sig

        # plt.subplot(1,3,3)
        # plt.plot(self.working_audio)
        # plt.show()

        sd.play(self.working_audio, self.fs)
        #sd.wait(self.dur+0.5)


    def convert_to_original(self):
        """ Convert back to original audio data type """
        # 1. Multiply float64 by original data type max
        sig = self.working_audio * self.wav_dict[str(self.data_type)][1]
        if self.data_type != 'float32':
            # 2. Round to return to integer values
            sig = np.round(sig)
        # 3. Convert back to original data type
        sig = sig.astype(self.data_type)
        print(f"Converted data type: {str(type(sig[0]))}")
        self.working_audio = sig


    @staticmethod
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


    @staticmethod
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


    def rms(self):
        """ 
            Calculate the root mean square of a signal. 
            
            NOTE: np.square will return invalid, negative 
                results if the number excedes the bit 
                depth. In these cases, convert to int64
                EXAMPLE: sig = np.array(sig,dtype=int)

            Written by: Travis M. Moore
            Last edited: Feb. 3, 2020
        """
        theRMS = np.sqrt(np.mean(np.square(self.working_audio)))
        return theRMS


    def setRMS(self, eq='n'):
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
        amp = self.level
        sig = self.working_audio
        if len(sig.shape) == 1:
            rmsdb = self.mag2db(self.rms())
            refdb = amp
            diffdb = np.abs(rmsdb - refdb)
            if rmsdb > refdb:
                sigAdj = sig / self.db2mag(diffdb)
            elif rmsdb < refdb:
                sigAdj = sig * self.db2mag(diffdb)
            # Edit 5/17/22
            # Added handling for when rmsdb == refdb
            elif rmsdb == refdb:
                sigAdj = sig
            return sigAdj
            
        elif len(sig.shape) == 2:
            rmsdbLeft = self.mag2db(self.rms(sig[0]))
            rmsdbRight = self.mag2db(self.rms(sig[1]))

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
                sigAdjLeft = sig[0] / self.db2mag(diffdbLeft)
            elif rmsdbLeft < refdb:
                sigAdjLeft = sig[0] * self.db2mag(diffdbLeft)
            # Adjust right channel
            if rmsdbRight > refdb:
                sigAdjRight = sig[1] / self.db2mag(diffdbRight)
            elif rmsdbRight < refdb:
                sigAdjRight = sig[1] * self.db2mag(diffdbRight)

            # If there is a lvl difference to maintain across channels
            if eq == 'n':
                if lvlAdv == 'left':
                    sigAdjLeft = sigAdjLeft * self.db2mag(ILD/2)
                    sigAdjRight = sigAdjRight / self.db2mag(ILD/2)
                elif lvlAdv == 'right':
                    sigAdjLeft = sigAdjLeft / self.db2mag(ILD/2)
                    sigAdjRight = sigAdjRight * self.db2mag(ILD/2)

            sigBothAdj = np.array([sigAdjLeft, sigAdjRight])
            return sigBothAdj
