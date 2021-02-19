# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 13:32:33 2018

Simple script to copy and rename files from hyperscanning backup onto my laptop

@author: adamb
"""


#%% IMPORTS

import argparse
import os
import wave
import glob
from shutil import copy2


#%% Wave creator function for the binary recorded data

def wavmaker(filename, CHANNELS, RATE):

    # create wav for audio
    f = open(filename, 'rb')
    audio = f.read()
    wavef = wave.open(filename+'.wav', 'w')
    wavef.setnchannels(CHANNELS)
    wavef.setsampwidth(2)
    wavef.setframerate(RATE)
    wavef.writeframes(audio)
    wavef.close()
    f.close()

    return


#%% File copy, etc

def main(pairN):

    # source path
    srcPath = '/media/adamb/TOSHIBA_EXT/hyperscanning_data_backup/pair'\
              + str(pairN)
    # target path
    tgPath = '/home/adamb/Documents/hyperscanning_pairData/pair'\
             + str(pairN)

    # if target path doesn't exist yet, create folder
    if not os.path.exists(tgPath):
        os.mkdir(tgPath)

    # same for the DBIC and DHMC files
    for site in ['DBIC', 'DHMC']:

        # check if we need account for an extra folder layer
        if not glob.glob(os.path.join(srcPath, site, 'behav/run1*')):
            extra = '*'
        else:
            extra = ''

        # go through both runs
        for run in ['run1*', 'run2*']:

            # check if run1 is joint or individual
            if 'JointStory' in glob.glob(os.path.join(srcPath, site,
                                                      'behav', extra,
                                                      run))[0]:
                story = 'JS'
            else:
                story = 'IS'

            # if wav does not exist, create it
            if not glob.glob(os.path.join(srcPath, site, 'behav', extra,
                                          run, 'RecordedAudio.wav')):
                if glob.glob(os.path.join(srcPath, site, 'behav', extra,
                                          run, 'RecordedAudio')):
                    file = glob.glob(os.path.join(srcPath,
                                                  site,
                                                  'behav',
                                                  extra,
                                                  run,
                                                  'RecordedAudio'))[0]
                    wavmaker(file, 1, 16000)
                else:
                    print('Missing audio file? Cannot ' +
                          'find the \'RecordedAudio\' file for pair ' +
                          str(pairN) + ', for ' + site)

            # copy and rename wav file
            if run == 'run1*':
                tgFile = tgPath + '/RecordedAudio_' + site +\
                         '_pair' + str(pairN) + '_cond1_' + story + '.wav'
            else:
                tgFile = tgPath + '/RecordedAudio_' + site +\
                         '_pair' + str(pairN) + '_cond2_' + story + '.wav'
            copy2(glob.glob(os.path.join(srcPath, site, 'behav', extra, run,
                                         'RecordedAudio.wav'))[0],
                  tgFile)


#%% MAIN

if __name__ == '__main__':

    # input arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'pairNumber',
        type=int,
        default=0,
        help='pair number (1-14)')
    # parse
    args = parser.parse_args()
    main(args.pairNumber)
