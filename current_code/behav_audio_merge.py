"""
FOR DBIC-DHMC DATA

USAGE: python3 behav_audio_merge -i BEHAV_DATA_FOLDER


See parameters / magic numbers in the first function "params".

The script does two things:

(1) We first look up all audio files that only exist in binary form, 
not in wav files yet. Specifically, we glob.glob() recursively
in the input folder for files named "RecordedAudio", "PlayedAudio",
and "ReceivedAudio". If a file does not have a corresponding wav file,
a wav file is created for it with the same path and filename, just adding
the ".wav" extension to it. This part is handled by the 
"wav_from_binary_audio" and "wavmaker" functions.

(2) We then create stereo wav files for each pair and run specified
in the "params" function, by combining the two "RecordedAudio.wav" files
as two mono channels into one file. The resulting wav file is named
"pair(PAIRNUMBER)_run(RUNNUMBER)_stereo.wav" and is saved into 
the input_folder.
For combining the two mono channels, we also align them in time. This
is achieved by looking up the audio file byte positions that 
corresponded to the common start time of the given pair and run. 
This information was saved out into the "TimingsLog.txt" files -
see the fifth line starting with "fOut:". The number there was the file
position of the "RecordedAudio" binary at the common start time. 
This is transformed into frame number by simply dividing by two, as we
used paInt16 format for audio that coded each frame into 16 bits = 2 bytes.
Using these byte / frame positions, the "RecordedAudio.wav" channels are 
trimmed to their start positions at common start time, and then combined
into the stereo wav file.
This part is handled by the functions "find_wavs_bytepos" and 
"get_Stereo_wavs".
 

"""

# Imports
from os import path
import argparse
import glob
import wave
import sys
import numpy as np
from scipy.io import wavfile


# magic numbers, parameters
def params():
    search_terms = ['RecordedAudio', 'ReceivedAudio', 'PlayedAudio']
    audio_rate = 16000
    pair_range = [1, 14]
    run_range = [1, 2]
    stereo_wav_term = 'RecordedAudio'
    log_file_term = 'TimingsLog.txt'
    return search_terms, audio_rate, \
           pair_range, run_range, stereo_wav_term, \
           log_file_term


# Wave creator function for the binary recorded data
def wavmaker(filename, channels, rate):
    # create wav for audio
    f = open(filename, 'rb')
    audio = f.read()
    wavef = wave.open(filename+'.wav', 'w')
    wavef.setnchannels(channels)
    wavef.setsampwidth(2)
    wavef.setframerate(rate)
    wavef.writeframes(audio)
    wavef.close()
    f.close()
    return


# find binary audio files matching the search_terms and create wav files from them
def wav_from_binary_audio(folder, search_terms, audio_rate):
    # go through all given glob prompts
    for term in search_terms:
        tmp_list = glob.glob(folder + '**/*' + term, recursive=True)
        # only include files in the results list which have no corresponding wav file already
        binary_audio_list = [f for f in tmp_list if not path.exists(f + '.wav')]
        # create wav from binary audio files
        for binary_audio in binary_audio_list:
            wavmaker(binary_audio, channels=1, rate=audio_rate)
    return


# find the corresponding ReceivedAudio wav files for a given pair and run
def find_wavs_bytepos(folder, pair_no, run_no, wav_term, timestamp_term):

    # Look for the wav files with glob
    print('\nLooking for ' + wav_term + '.wav files for pair ' + str(pair_no) +
          ', run ' + str(run_no) + '...')
    tmp_search_term = folder + '**/pair' + str(pair_no) + '/**/run' + str(run_no) + '*/' + wav_term + '.wav'
    wavs = glob.glob(tmp_search_term, recursive=True)
    # if there is a problem with the number of wav files, quit
    if len(wavs) != 2:
        print('Found too many or not enough wav files!')
        print(wavs)
        sys.exit()
    else:
        print('Found wav files: ')
        print(wavs)

    # Look for the timestamp log files
    print('\nLooking for log files...')
    tmp_search_term = folder + '**/pair' + str(pair_no) + '/**/run' + str(run_no) + '*/' + timestamp_term
    log_files = glob.glob(tmp_search_term, recursive=True)
    # if there is a problem with the number of log files, quit
    if len(log_files) != 2:
        print('Found too many or not enough log files!')
        print(log_files)
        sys.exit()
    else:
        print('Found log files: ')
        print(log_files)

    # extract byte position of audio files at common start
    bytepos = []
    start_times = []
    for i in range(2):
        with open(log_files[i], 'r') as fopen:
            txt = fopen.readlines()
        if txt[4][0:5] == 'fOut:':
            bytepos.append(txt[4][6:].strip())
            start_times.append(txt[0][17:].strip())
    print('\nExtracted common start times from logs:' +
          '\n' + start_times[0] +
          '\n' + start_times[1])
    print('\nExtracted audio file positions (in bytes) at common start time:' +
          '\n' + bytepos[0] +
          '\n' + bytepos[1])

    bytepos = [int(f) for f in bytepos]

    return wavs, bytepos


# # Function to trim two wav files according to the received starting byte positions,
# # and also to the same length
# def trim_wavs(wavs, bytepos, folder, pair_no, run_no):
#     # Read in wav files and set their file positions to the ones at common start time
#     print('\nReading in wave files')
#     print(wavs)
#     channel1_wave = wave.open(wavs[0])
#     channel2_wave = wave.open(wavs[1])
#     channel1_wave.setpos(bytepos[0])
#     channel2_wave.setpos(bytepos[1])
#     print('\nWave file positions set to common start time positions')
#
#     # Open a wav file for writing
#     stereo_filename = folder + 'pair' + str(pair_no) + '_run' + str(run_no) + '_stereo.wav'
#     stereo_wav = wave.open(stereo_filename,'w')
#
#     # set the parameters of stereo wav
#     nchannels = 2
#     sampwidth = channel1_wave.getsampwidth()
#     framerate = channel1_wave.getframerate()
#     comptype = channel1_wave.getcomptype()
#     compname = channel1_wave.getcompname()
#     # The length of the stereo audio is defined by
#     # the length of the mono audios,
#     # minus the effect of trimming in the beginning
#     nframes = min(channel1_wave.getnframes() - bytepos[0] / 2, channel2_wave.getnframes() - bytepos[1] / 2)
#     stereo_wav.setparams((nchannels, sampwidth, framerate, nframes,
#                         comptype, compname))
#
#     # write frame into stereo wav:
#     # in the case of stereo, interlaced frames are required
#
#     return


def get_stereo_wavs(wavs, bytepos, folder, pair_no, run_no):
    # Read in wav files, trim them
    print('\nReading in wave files')
    print(wavs)
    frate1, channel1 = wavfile.read(wavs[0])
    frate2, channel2 = wavfile.read(wavs[1])
    if frate1 != frate2:
        print('\nWav files have different frame rates!')
        sys.exit()
    # byte position to frames, then trim from the beginning
    channel1 = channel1[int(bytepos[0]/2):]
    channel2 = channel2[int(bytepos[1]/2):]
    # get min remaining size and trim again
    minsize = min(channel1.size, channel2.size)
    channel1 = channel1[:minsize]
    channel2 = channel2[:minsize]
    # create matrix for stereo output wav
    channels = np.vstack([channel1, channel2])

    # create stereo output wav
    stereo_filename = folder + 'pair' + str(pair_no) + '_run' + str(run_no) + '_stereo.wav'
    print('\nWriting stereo wav file to ' + stereo_filename)
    print(channels.shape)
    wavfile.write(stereo_filename, frate1, channels.T)

    return



def main():

    # parse input arguments
    parser = argparse.ArgumentParser()

    # Input arg "input_folder"
    parser.add_argument(
        '-i',
        '--input_folder',
        type=str,
        help='Folder of behavioral (audio) data.'
             'Glob is used to find all files that need to be transformed to wav.'
             'Glob is recursive!')

    # parse arguments, get list
    args = parser.parse_args()

    # check inputs
    if not path.exists(args.input_folder):
        raise ValueError('Input arg "input_folder" is not found')

    # start message
    print('\nbehav_audio_merge was started with input_folder: ',
          args.input_folder)

    # call params
    search_terms, audio_rate, pair_range, \
        run_range, stereo_wav_term, \
        log_file_term = params()

    # get list of files we need to transform to wav
    wav_from_binary_audio(args.input_folder, search_terms, audio_rate)

    # for each pair and run, find the two corresponding RecordedAudio wav files,
    # and the corresponding audio frame numbers at common start time
    for pair_no in range(pair_range[0], pair_range[1]+1):
        for run_no in range(run_range[0], run_range[1]+1):
            wavs, bytepos = find_wavs_bytepos(args.input_folder,
                                              pair_no,
                                              run_no,
                                              stereo_wav_term,
                                              log_file_term)
            # load up the two wav files and trim them so that
            # they start from the common start time and end
            # at the same time too
            get_stereo_wavs(wavs, bytepos, args.input_folder, pair_no, run_no)

    return


#
if __name__ == '__main__':
    main()
