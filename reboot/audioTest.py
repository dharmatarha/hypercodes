#!/usr/bin/env python3.9
# -*- coding: utf-8 -*-
"""

@author: adamb

Minimal audio-only version.

"""

# imports
import argparse
import socket
import pyaudio
import sys
import time
import multiprocessing


#%% MAGIC NUMBERS aka hardcoded variables, default values
def magicNumbers():
    # default audio settings: mono, 16kHz sampling, 16bit encoding
    CHANNELS = 1
    RATE = 16000
    FORMAT = pyaudio.paInt16
    # default filenames (for saving audio)
    savefileOut = 'RecordedAudio'
    savefileIn1 = 'ReceivedAudio'
    savefileIn2 = 'PlayedAudio'
    savefileLog = 'TimingsLog.txt'
    savefileTTL = 'TTLtimestamps.txt'
    # default port numbers
    # local
    portIn = 30302
    portOut = 30301
    # remote
    PortIn = 30302
    PortOut = 30301
    # keys to check at various times
    keylist = ['escape']

    return(CHANNELS, RATE, FORMAT, savefileOut,
           savefileIn1, savefileIn2, savefileLog, savefileTTL, portIn, portOut,
           PortIn, PortOut, keylist)


#%% Socket function
# Opens simple UDP socket, binds it to given port number at localhost.

def openSocket(port):
    socketFlag = False
    # define socket
    try:
        socketUDP = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        print('\nSocket created')
    except socket.error:
        print('\nFailed to create UDP socket')
    # bind port
    host = ''  # localhost
    try:
        socketUDP.bind((host, port))
        print('\nUDP socket bound to local port: ', port)
        socketFlag = True
    except socket.error:
        print('\nFailed to bind UDP socket to local port ', port)
    return socketFlag, socketUDP


#%% Callback function for non-blocking pyaudio (portaudio) input. Minimal
# version.
def callbackInput(in_data, frame_count, time_info, status):
    # send new data to other side
    try:
        socketOut.sendto(in_data, (IP, PortIn))
    except socket.error:
        print('Could not send a packet\n')
    # return data and flag
    return (in_data, pyaudio.paContinue)


#%% Callback function for non-blocking pyaudio (portaudio) output.

# In this version, we use a simple continuous buffer that collects
# all incoming packets on one end and is read out by the callback on the other.
# Minimal version for testing.
def callbackOutput(in_data, frame_count, time_info, status):
    global lastData
    # once the buffer is filled for the first time, startFlag is set and
    # callback can read from it
    if startFlag:
        # first check if there is enough data available to read
        if len(audioBuffer) > CHUNK*2:
            data = audioBuffer[0:CHUNK*2]
            del audioBuffer[0:CHUNK*2]
            lastData = data
        # if buffer is empty, update the underflow counter
        else:
            data = lastData

    # until startFlag is set, callback reads from a silence buffer (zeros)
    else:
        if len(silenceBuffer) > CHUNK*2:
            data = silenceBuffer[0:CHUNK*2]
            del silenceBuffer[0:CHUNK*2]
            lastData = data
        else:
            data = lastData

    data = bytes(data)
    return data, pyaudio.paContinue


#%% Function to set up mic
# Uses pyaudio (portaudio) for a non-blocking input device.
# Device is default input set on platform.

def micOpen(FORMAT, CHANNELS, RATE, CHUNK):
    p = pyaudio.PyAudio()
    stream = p.open(format=FORMAT,
                    channels=CHANNELS,
                    rate=RATE,
                    input=True,
                    frames_per_buffer=CHUNK,
                    stream_callback=callbackInput,
                    start=False)  # IMPORTANT: don't start yet
    return stream, p


#%% Function to open output device
# Uses pyaudio (portaudio). Chooses default output device on platform.

def speakersOpen(FORMAT, CHANNELS, RATE, CHUNK):
    # open pyaudio (portaudio) device
    p = pyaudio.PyAudio()
    # open portaudio output stream
    stream = p.open(format=FORMAT,
                    channels=CHANNELS,
                    rate=RATE,
                    output=True,
                    frames_per_buffer=CHUNK,
                    stream_callback=callbackOutput,
                    start=False)  # IMPORTANT: don't start yet
    return stream, p


#%% Main input (microphone) function. Handles audio input and
# transmission. Should be called in separate process (multiprocessing!),

def inputProcess(FORMAT, CHANNELS, RATE, CHUNK, queueInput):
    global chunkCounter, streamInput, pIn
    # open input dev
    streamInput, pIn = micOpen(FORMAT, CHANNELS, RATE, CHUNK)
    # print start message
    print('\nEverything seems all right, input stream open on our side.\n')

    # start input stream
    streamInput.start_stream()

    # kepp channel open until asked to quit
    while streamInput.is_active():
        time.sleep(0.01)
        # if escape key was pressed, terminate
        if not queueInput.empty():
            break

    # input cleanup
    print('\nTransmission finished, cleaning up input...\n')
    # pyaudio
    streamInput.stop_stream()
    streamInput.close()
    pIn.terminate()
    # sockets
    socketOut.close()
    print('Closed input stream, socket\n')

    return


#%% Main output (receiver) function. Handles audio output and
# packet control. Should be called in separate process (multiprocessing!),
# after networkInit(), at the same time as inputProcess()

def outputProcess(BUFFER, CHUNK, FORMAT, CHANNELS,
                  RATE, queueOutput):

    # these need to be global...
    global startFlag, audioBuffer, streamOutput
    global pOut, silenceBuffer, lastData

    # initialize callback start flag
    startFlag = 0

    # counter for all received UDP packets
    packetCounter = 0

    # open output dev
    streamOutput, pOut = speakersOpen(FORMAT, CHANNELS, RATE, CHUNK)
    print('\nAudio output set up, waiting for transmission.')

    # create buffer for incoming packets
    audioBuffer = bytearray()

    # start stream with a silent buffer (silence buffer)
    silenceBuffer = b'x\00'*2*CHUNK*BUFFER
    silenceBuffer = bytearray(silenceBuffer)
    lastData = silenceBuffer[0:CHUNK*2]
    streamOutput.start_stream()

    # Main loop: listen for UDP, fill buffer, hand it to output stream
    while True:
        # if escape key was pressed, ,terminate
        if not queueOutput.empty():
            break
        # receive UDP packet - remember this is in non-blocking mode now!
        packet = []
        try:
            packet = socketIn.recv(CHUNK*4)
        except:
            pass
        # if we received anything
        if packet:
            # adjust packet counter
            packetCounter += 1
            # put data into buffer
            audioBuffer.extend(packet)

            # set startFlag for callback once buffer is filled for the first
            # time
            if packetCounter == BUFFER:
                startFlag = 1

            # if audioBuffer is getting way too long, chop it back, the
            # treshold is two times the normal size
            if len(audioBuffer) > 2*CHUNK*BUFFER*3:
                del audioBuffer[0:2*CHUNK*BUFFER]

            # display state
            if packetCounter % 250 == 0:
                print('Chunk no ', packetCounter, 'received successfully')
                print('Current buffer size: '+str(len(audioBuffer)))

    # cleanup
    print('\nTransmission finished, cleaning up output...\n')
    # pyaudio
    streamOutput.stop_stream()
    streamOutput.close()
    pOut.terminate()
    # sockets
    socketIn.close()
    print('Output stream socket closed\n')

    return


#%%
def audioAll():
    # sockets need to be global for callbacks
    global socketIn, socketOut, PortIn, PortOut
    # call hardcoded variables
    [CHANNELS, RATE, FORMAT, savefileOut,
     savefileIn1, savefileIn2, savefileLog, savefileTTL, portIn, portOut,
     PortIn, PortOut, keylist] = magicNumbers()
    # UPD sockets for transmission
    socketFlag, socketOut = openSocket(portOut)
    if not socketFlag:
        print('\n\nCould not create or bind UDP socket. Wtf.')
        sys.exit()
    socketFlag, socketIn = openSocket(portIn)
    if not socketFlag:
        print('\n\nCould not create or bind UDP socket. Wtf.')
        sys.exit()
    # set sockets to non-blocking
    socketOut.settimeout(0)
    socketIn.settimeout(0.1)

    # audio I/O processes run in separate processes
    queueInput = multiprocessing.Queue()
    queueOutput = multiprocessing.Queue()
    audioInput = multiprocessing.Process(name='audioInput',
                                         target=inputProcess,
                                         args=(FORMAT,
                                               CHANNELS,
                                               RATE,
                                               CHUNK,
                                               queueInput,))
    audioOutput = multiprocessing.Process(name='audioOutput',
                                          target=outputProcess,
                                          args=(BUFFER,
                                                CHUNK,
                                                FORMAT,
                                                CHANNELS,
                                                RATE,
                                                queueOutput,))
    audioInput.start()
    audioOutput.start()
    print('Audio devices started in separate proesses\n')

    # keep channels open for some time
    start = time.time()
    while time.time()-start < 30:
        time.sleep(0.5)

    # close processes, clean up
    queueOutput.put('die')
    queueInput.put('die')
    audioOutput.join()
    audioInput.join()


#%% Main
if __name__ == '__main__':

    # input arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-c',
        '--CHUNK',
        nargs='?',
        type=int,
        default=512,
        help='Audio chunk (packet) size in frames (1 frame = 2 bytes with ' +
             'current format settings). Integer. Default = 512')
    parser.add_argument(
        '-i',
        '--IP',
        nargs='?',
        type=str,
        default='',
        help='IP (ipv4) in string for local network transmission. ' +
             'Default = empty string, that is, localhost.')
    parser.add_argument(
        '-b',
        '--BUFFER',
        nargs='?',
        type=int,
        default=4,
        help='No. of chunks to buffer for audio output. Integer. Default = 4')
    args = parser.parse_args()

    # check inputs
    # the following check for power of two is a really cool trick!
    if ((args.CHUNK != 0) and ((args.CHUNK & (args.CHUNK - 1)) == 0) and
       (not (args.CHUNK < 128)) and (not (args.CHUNK > 4096))):
        pass
    else:
        raise ValueError('CHUNK should be power of two, between 128 and 4096.')
    if 1 <= args.BUFFER <= 25:
        pass
    else:
        raise ValueError('Unexpected value for argument BUFFER. ',
                         '(please have it 1 <= and <= 25.')

    # some global -local handling:
    # anything used later in callbacks needs to be global
    global IP, BUFFER, CHUNK
    IP = args.IP
    BUFFER = args.BUFFER
    CHUNK = args.CHUNK

    audioAll()

    # End
    print('\n\nEverything ended / closed the way we expected. Goodbye!\n\n')
