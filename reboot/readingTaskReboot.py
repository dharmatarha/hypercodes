#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Simple reading task. The display is refreshed with new text periodically,
participants are asked to simply read what they see.

Works keypress style TTLs (simple key presses of "5" via an usb keyboard-like input device)

Input args:.
--pair (-p): Pair number, integer, between 1-99. Defaults to 99. 
             It is added to the log file name.
--site (-s): Site name, string, one of ["DBIC", "Harvard", "Test"]. Defaults to "Test". 
             It is added to the log file name.

Outputs a log file of all major things that happen, this includes DHMC-style
TTLs.

Created Nov 29 2022

@author: adamb
"""

#%% Imports

import wave
import time
import pyaudio
from io import open
from psychopy import visual, core, event
import argparse
import multiprocessing
import serial
import codecs
import sys


#%% magic numbers, parameters
def magicNumbers():
    # audio chunk
    audioChunk = 1024
    # default audio settings: mono, 16kHz sampling, 16bit encoding
    audioChannels = 1
    audioRate = 16000
    audioFormat = pyaudio.paInt16
    # number of TTL pulses to wait before starting audio recording
    nTTL = 3
    # read in the stimulus we will use
    f = codecs.open('stimFinal.txt', 'r', encoding='utf8')
    text = f.readlines()
    f.close()
    # read in the corresponding syllable lengths, one integer per text line
    f = codecs.open('sylFinal.txt', 'r', encoding='utf8')
    syls = f.readlines()
    sylN = [int(syl) for syl in syls]
    f.close()
    # filenames for recorded audio and logs
    fOutName = 'recordedSpeech'
    logName = 'readingLog'
    # time for one syllable, in frames (we assume 60 Hz and 333 ms / syl)
    sylFrames = 20
    # keys to check at various times
    keylist = ['5', 'escape']

    return(audioChunk, audioChannels, audioRate, audioFormat, nTTL,
           text, sylN, fOutName, logName, sylFrames, keylist)


#%% Function to set up mic
# Uses pyaudio (portaudio) for a non-blocking input device.
# Device is default input set on platform.

def micOpen(audioFormat, audioChannels, audioRate, audioChunk):
    p = pyaudio.PyAudio()
    stream = p.open(format=audioFormat,
                    channels=audioChannels,
                    rate=audioRate,
                    input=True,
                    frames_per_buffer=audioChunk,
                    stream_callback=callbackInput,
                    start=False)  # IMPORTANT: don't start yet
    return stream, p


#%% Callback function for non-blocking pyaudio (portaudio) input
# Callback saves input audio and handles timestamps too.
# IMPORTANT:
# Expects output file to be open and ready for writing.

def callbackInput(in_data, frame_count, time_info, status):

    # keep track of chunks
    global chunkCounter
    # refresh counter, keep track of progress
    chunkCounter += 1
    if chunkCounter % 500 == 0:
        print('\nAudio chunk # ' + str(chunkCounter) + ' recorded\n')
    # write out new data before we mess with it
    fOut.write(in_data)
    # return data and flag
    return (in_data, pyaudio.paContinue)


#%% Main input (microphone) function. Handles audio input and
# transmission. Should be called in separate process (multiprocessing). Relies
# on corresponding callback function

def inputProcess(audioFormat, audioChannels,
                 audioRate, audioChunk, queueInput):

    global chunkCounter, fOut
    # init chunkCounter
    chunkCounter = 0

    # open input dev
    streamInput, pIn = micOpen(audioFormat,
                               audioChannels,
                               audioRate,
                               audioChunk)

    # print start message
    print('\nEverything seems all right, mic is listening.\n')

    # start input stream
    start = time.time()
    streamInput.start_stream()

    # wait until all audio is sent
    while streamInput.is_active():
        time.sleep(0.01)
        # if escape key was pressed, terminate
        if queueInput.get() == 'die':
            break

    # message
    print('\n' + str(chunkCounter)+' chunks saved, time taken: ' +
          str(time.time()-start) + '\n')

    # input cleanup
    streamInput.stop_stream()
    streamInput.close()
    pIn.terminate()
    # files
    fOut.close()

    return


#%% Wave creator function for the binary recorded data
def wavmaker(filename, audioChannels, audioRate):
    # create wav for audio
    f = open(filename, 'rb')
    audio = f.read()
    wavef = wave.open(filename+'.wav', 'w')
    wavef.setnchannels(audioChannels)
    wavef.setsampwidth(2)
    wavef.setframerate(audioRate)
    wavef.writeframes(audio)
    wavef.close()
    f.close()
    return


#%% Experimental function integrating everything

def GO(pair, site):

    # declare the recording file as global - we will use it in the callback
    global fOut

    # load parameters
    [audioChunk, audioChannels, audioRate, audioFormat, nTTL,
     text, sylN, fOutName, logName, sylFrames, keylist] = magicNumbers()
    # add pair number and site name to output file names
    fOutName = fOutName + '_pair' + str(pair) + '_' + site
    logName = logName+ '_pair' + str(pair) + '_' + site + '.txt'
    # init audio file
    fOut = open(fOutName, 'wb')
    # init logfile
    GO_reading = open(logName, 'w')
    print('Log file opened, ' + format(time.time(), '.4f'),
          file=GO_reading)

    # start clock for recording key presses
    startClock = core.Clock()
    # get the timestamp for starting the clock
    startStamp = startClock.getLastResetTime()
    print('\nStartClock reset time, ' + format(startStamp, '.4f'),
          file=GO_reading)

    # set up window: black, fullscreen, try a common size for projectors
    try:
        win = visual.Window([1600, 900],
                            color='black',
                            monitor='testMonitor',
                            fullscr=True)
        print('\nMonitor frame period: ' + str(win.monitorFramePeriod)) 
    except:
        print('\nProblem while setting up window\n')
        sys.exit()

    # audio input runs in separate process
    queueInput = multiprocessing.Queue()
    audioInput = multiprocessing.Process(name='audioInput',
                                         target=inputProcess,
                                         args=(audioFormat,
                                               audioChannels,
                                               audioRate,
                                               audioChunk,
                                               queueInput,))
    audioInput.start()
    print('Audio input process started. We are recording, ' +
          format(time.time(), '.4f'), file=GO_reading)

    # Instructions screen, with start message
    instText = ('Read out loud the two short stories you will be presented '
                'with, line by line, as if you were telling it yourself.\n\n'
                'A moving dot will indicate the approximate speed to '
                'read with.\n\n'
                'Try not to move your head too much '
		'and please use a normal speaking voice, '
		'even if it is hard to hear yourself in the scanner.\n\n'
                )
    instruction = visual.TextStim(win,
                                  instText,
                                  color='white',
                                  height=0.08)
    instruction.draw()
    # flip instruction onto screen, get timestamp
    win.flip()
    instFlip = time.time()
    print('\nFlipping instructions, ' +
          format(instFlip, '.4f'), file=GO_reading)

    # Instructions amended with the Starting... part
    instText = ('Read out loud the two short stories you will be presented '
                'with, line by line, as if you were telling it yourself.\n\n'
                'A moving dot will indicate the approximate speed to '
                'read with.\n\n'
                'Try not to move your head too much '
		'and please use a normal speaking voice, '
		'even if it is hard to hear yourself in the scanner.\n\n'
                'Waiting for the scanner, just a moment...'
                )
    instruction = visual.TextStim(win,
                                  instText,
                                  color='white',
                                  height=0.08,
                                  alignVert='center',
                                  alignHoriz='center')
    instruction.draw()

    # wait a bit before switching, listen to keys in the meantime
    while time.time() < (instFlip + 10):
        core.wait(0.01)
        # capture key presses, timing is relative to keyClock
        keys = event.getKeys(keylist,
                             timeStamped=startClock)
        if keys:
            # if event.getKeys returns an escape, we break the while loop
            if keys[0][0] == 'escape':
                print('\nEscape was pressed, ' +
                      format(keys[0][1] + startStamp, '.4f')+'\n',
                      file=GO_reading)
                break
    win.flip()

    # Wait for three ttl-s in the form of keypresses
    counter = 0
    while counter < 3:
        core.wait(0.01)
        # capture key presses, timing is relative to keyClock
        keys = event.getKeys(keylist,
                             timeStamped=startClock)
        if keys:
            if keys[0][0] == '5':
                print('\nTTL, ' +
                      format(keys[0][1] + startStamp, '.4f')+'\n',
                      file=GO_reading)
                counter += 1
            # if event.getKeys returns an escape, we break the while loop
            elif keys[0][0] == 'escape':
                print('\nEscape was pressed, ' +
                      format(keys[0][1] + startStamp, '.4f')+'\n',
                      file=GO_reading)
                break

    # START!
    win.flip()
    print('\nTTL Start signal arrived. Starting main loop, ' +
          format(time.time(), '.4f'), file=GO_reading)

    # prepare main loop
    emptyFlag = 0
    escFlag = 0
    counter = 0
    winW, winH = win.size
    # Loop through lines to read out loudly
    for line in text:
        if escFlag:
            break
        # text stimulus
        stim = visual.TextStim(
                               win,
                               line,
                               color='white',
                               height=0.14
                               )
        # get the dimensions of the text box
        sW, sH = stim.boundingBox
        # the bounding box is given in pixels, we calculate it in normalized
        sWRatio = sW/winW
        sHRatio = sH/winH
        # calculate the width value for repositioning the rect on each frame
        wUpdate = 2*sWRatio / (sylFrames*sylN[counter])
        # dot stimulus
        rect = visual.Rect(
                           win=win,
                           width=0.02,
                           height=0.03,
                           fillColor='gray',
                           lineColor='gray',
                           pos=[-2*sWRatio/2, 2*sHRatio+0.02]
                           )

        # if empty line
        if line == '\n':
            emptyFlag = 1

        for frameN in range(sylFrames*sylN[counter]):
            rect.pos += (wUpdate, 0)
            stim.draw()
            if not emptyFlag:
                rect.draw()
            win.flip()
            # capture key presses, timing is relative to keyClock
            keys = event.getKeys(keylist,
                                 timeStamped=startClock)
            if keys:
                if keys[0][0] == '5':
                    print('\nTTL, ' +
                          format(keys[0][1] + startStamp, '.4f')+'\n',
                          file=GO_reading)
                # if event.getKeys returns an escape, we break the while loop
                elif keys[0][0] == 'escape':
                    print('\nEscape was pressed, ' +
                          format(keys[0][1] + startStamp, '.4f')+'\n',
                          file=GO_reading)
                    escFlag = 1
                    break
            # log timing at the first flip
            if frameN == 0:
                print('\nLine number ' + str(counter) + ' displayed, ' +
                      format(time.time(), '.4f'), file=GO_reading)
        counter += 1
        emptyFlag = 0

    # last line display ended
    print('Last line is over ' +
          format(time.time(), '.4f'), file=GO_reading)

    # show end message
    endText = ('The task has ended, great job!')
    endMess = visual.TextStim(win,
                              endText,
                              color='white',
                              height=0.14)
    endMess.draw()
    win.flip()
    print('\nEnd message displayed ' +
          format(time.time(), '.4f'), file=GO_reading)
    core.wait(2)

    # close everything
    queueInput.put('die')
    audioInput.join()
    print('\nAudio input process ended and joined, ' +
          format(time.time(), '.4f'), file=GO_reading)
    print('\nAudio input process ended.\n')

    # close log file
    GO_reading.close()
    print('\nClosed log file.\n')

    # create wave from recorded audio
    wavmaker(fOutName, audioChannels, audioRate)
    print('\nSaved wav file.\n')

    return


#%% Main

if __name__ == '__main__':

    # input arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-p',
        '--pair',
        nargs='?',
        type=int,
        default=99,
        help='Pair number, integer. Defaults to 99.')
    parser.add_argument(
        '-s',
        '--site',
        nargs='?',
        type=str,
        default='Test',
        help='Scanner site, string. "DBIC", "Harvard", or "Test". Defaults to "Test".')
    args = parser.parse_args()

    # check input
    if 1 <= args.pair <= 99:
        pass
    else:
        raise ValueError('Unexpected value for argument --pair ',
                         '(should be in range 1 - 99)')
    if args.site in ['DBIC', 'Harvard', 'Test']:
        pass
    else:
        raise ValueError('Unexpected value for argument --site ',
                         '(should be "DBIC" or "Harvard")')


    GO(args.pair, args.site)
