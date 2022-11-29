#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simple audio script.
Plays one wave file while displaying instructions, records TTLs.

Works with keypress style TTLs (simple key presses of "5" via an usb keyboard-like input device)

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
import sys


#%% magic numbers, parameters
def magicNumbers():
    # audio chunk
    chunk = 1024
    # read file
    audio = wave.open('audioStimNew.wav', 'rb')
    # number of TTL pulses to wait before starting audio output
    nTTL = 3
    # keys to check at various times
    keylist = ['5', 'escape']

    return(chunk, audio, nTTL, keylist)


#%% audio output stream
def audioStream(chunk, audio):
    # create device
    p = pyaudio.PyAudio()
    # create output stream, with params specified from audio object
    stream = p.open(
                    format=p.get_format_from_width(audio.getsampwidth()),
                    channels=audio.getnchannels(),
                    rate=audio.getframerate(),
                    output=True,
                    frames_per_buffer=chunk,
                    )
    # read in first chunk
    data = audio.readframes(chunk)

    return(p, stream, data)


#%% Cleanup
def cleanup(p, stream):
    stream.close()
    p.terminate()

    return


#%% Play audio
def GO(pair, site):

    # init logfile
    GO_audio = open('listeningLog_pair' + str(pair) + '_' + site + '.txt', 'w')
    print('\nScript started and log file opened, ' +
          format(time.time(), '.4f'),
          file=GO_audio)
    # load parameters
    chunk, audio, nTTL, keylist = magicNumbers()
    print('\nParameters loaded ' +
          format(time.time(), '.4f'),
          file=GO_audio)
    # init audio
    p, stream, data = audioStream(chunk, audio)
    print('\nAudio device started and first chunk prepared, ' +
          format(time.time(), '.4f'),
          file=GO_audio)

    # set up window: black, fullscreen, try a common size for projectors
    try:
        win = visual.Window([1366, 768],
                            color='black',
                            monitor='testMonitor',
                            fullscr=True)
    except:
        print('\nProblem while setting up window')
        sys.exit()

    # start clock for recording key presses
    startClock = core.Clock()
    # get the timestamp for starting the clock
    startStamp = startClock.getLastResetTime()
    print('\nStartClock reset time, ' + format(startStamp, '.4f'),
          file=GO_audio)

    # Instructions screen, with start message
    instText = ('Listen to the two short stories you are about to hear.\n\n'
                'During these stories, there will be a white dot displayed '
                'centrally. Please fixate on this dot while '
                'listening.\n\n'
                'We will not ask you to retell these stories, but please pay '
                'attention to them nonetheless\n\n'
                )
    instruction = visual.TextStim(win,
                                  instText,
                                  color='white',
                                  height=0.06)
    instruction.draw()
    # flip instruction onto screen, get timestamp
    win.flip()
    instFlip = time.time()

    # Instructions amended with the Starting... part
    instText = ('Listen to the two short stories you are about to hear.\n\n'
                'During these stories, there will be a white dot displayed '
                'centrally. Please fixate on this dot while '
                'listening.\n\n'
                'We will not ask you to retell these stories, but please pay '
                'attention to them nonetheless\n\n'
                'Waiting for the scanner, just a moment...'
                )
    instruction = visual.TextStim(win,
                                  instText,
                                  color='white',
                                  height=0.06,
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
                      file=GO_audio)
                break
    win.flip()

    # prepare simple fixation dot
    dotWin = visual.Circle(win,
                           units='pix',
                           radius=10,
                           lineColor='white')
    dotWin.fillColor = 'white'
    dotWin.draw()

    # Wait for three keypress type ttl-s (key "5", as at DHMC)
    counter = 0
    while counter < 3:
        core.wait(0.01)
        # capture key presses, timing is relative to startClock
        keys = event.getKeys(keylist,
                             timeStamped=startClock)
        if keys:
            if keys[0][0] == '5':
                print('TTL, ' +
                      format(keys[0][1] + startStamp, '.4f')+'\n',
                      file=GO_audio)
                counter += 1
            # if event.getKeys returns an escape, we break the while loop
            elif keys[0][0] == 'escape':
                print('\nEscape was pressed, ' +
                      format(keys[0][1] + startStamp, '.4f')+'\n',
                      file=GO_audio)
                break

    # switch to fixation dot, record time
    win.flip()
    print('\nFixation flip time:' + format(time.time(), '.4f'), file=GO_audio)

    # listen to keys while we wait a second
    dotWait = time.time()
    while time.time() < (dotWait + 1):
        core.wait(0.01)
        # capture key presses, timing is relative to startClock
        keys = event.getKeys(keylist,
                             timeStamped=startClock)
        if keys:
            # if event.getKeys returns an escape, we break the while loop
            if keys[0][0] == '5':
                print('TTL, ' +
                      format(keys[0][1] + startStamp, '.4f')+'\n',
                      file=GO_audio)
            elif keys[0][0] == 'escape':
                print('\nEscape was pressed, ' +
                      format(keys[0][1] + startStamp, '.4f')+'\n',
                      file=GO_audio)
                break

    # start time
    GO_start = format(time.time(), '.4f')
    print('\nAudio start time:' + GO_start, file=GO_audio)
    print('\n\nStarted audio at ' + GO_start)

    # play stream (looping from beginning of file to the end)
    while len(data) > 0:
        # writing to the stream is what *actually* plays the sound.
        stream.write(data)
        data = audio.readframes(chunk)
        # capture key presses, timing is relative to startClock
        keys = event.getKeys(keylist,
                             timeStamped=startClock)
        if keys:
            # if event.getKeys returns an escape, we break the while loop
            if keys[0][0] == '5':
                print('TTL, ' +
                      format(keys[0][1] + startStamp, '.4f')+'\n',
                      file=GO_audio)
            elif keys[0][0] == 'escape':
                print('\nEscape was pressed, ' +
                      format(keys[0][1] + startStamp, '.4f')+'\n',
                      file=GO_audio)
                break

    # wait till it stops
    stream.stop_stream()

    # end time
    GO_end = format(time.time(), '.4f')
    print('\nAudio end time:' + GO_end, file=GO_audio)
    print('\nLength in secs:' + str(float(GO_end)-float(GO_start)),
          file=GO_audio)
    print('\n\nAudio ended at ' + GO_end)
    print('Total time: ' + str(float(GO_end)-float(GO_start)))

    # show end message
    endText = ('The task has ended, great job!')
    endMess = visual.TextStim(win,
                              endText,
                              color='white',
                              height=0.06)
    endMess.draw()
    win.flip()
    print('\nEnd message displayed ' +
          format(time.time(), '.4f'), file=GO_audio)
    core.wait(2)


    # audio cleanup
    cleanup(p, stream)

    # psychopy cleanup
    core.wait(2)
    win.close()
    core.quit()

    # close log file
    GO_audio.close()
    print('\n\nLog file closed')


    print('\n\nCleaned up, closing')

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

