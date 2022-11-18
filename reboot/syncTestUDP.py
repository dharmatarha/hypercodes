#!/usr/bin/python3

import socket
import numpy as np
import time
import sys
import select
import argparse
import multiprocessing


def openSocket(port):
    """
    Opens an UDP port reachable on any address of the machine and binds it to specified port.
    Socket is set to non-blocking by default.

    Inputs
    port :          Numeric value, port number for binding socket

    Returns
    socketFlag :    Boolean, True if socket was created successfully
    socketUDP :     Socket obj, UDP, bound to input `port`, set to non-blocking
    """

    socketFlag = False
    socketCreated = False
    socketUDP = False
    # define socket
    try:
        socketUDP = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        socketUDP.setblocking(False)  # set to non-blocking
        socketCreated = True
        print('\nSocket created')
    except socket.error:
        print('\nFailed to create UDP socket')
    # bind port if it exists
    if socketCreated:
        host = ''  # "host address: '' represents INADDR_ANY, which is used to bind to all interfaces..."
        try:
            socketUDP.bind((host, port))
            print('\nUDP socket bound to port: ', port)
            socketFlag = True
        except socket.error:
            print('\nFailed to bind UDP socket to port ', port)
    return socketFlag, socketUDP


def closeSockets(socketList):
    """
    Helper closing sockets gracefully

    Inputs
    socketList :    List of socket objects
    """

    for s in socketList:
        s.close()

    return


def handshake(socketComm, addr):
    """
    Function to negotiate a shared starting time with other PC/process
    using an UDP socket for communication.

    Implements simple handshake procedure with two stages:
    (1) Send a predefined message to the supplied address ("addr") until there is an answer back
    (2) Send local time at start of stage two until there is incoming timestamp.
    Shared starting time is the average of the two time stamps plus a hardcoded constant ("commonStartDelay").

    Returns with "errorFlag" == True if
    (1) the socket is in error state at any point
    (2) any stage lasts longer than a hardcoded timeout constant ("handshakeTimeout")
    (3) local and remote timestamps differ more than a hardcoded constant ("timestampMaxDiff")

    Inputs
    socketComm :        Socket object, non-blocking UDP socket bound to a port
    addr :              Tuple containing remote address for UDP packets (IPv4, port)

    Returns
    errorFlag :         Boolean, True if handshake failed (see error conditions above)
    commonStartTime :   Float, negotiated starting time in Unix time. 0 if "errorFlag" == True
    """

    # default values for returns
    errorFlag = False
    commonStartTime = 0

    # flag for while loops
    successFlag = False

    # hardcoded constants
    handshakeTimeout = 60  # timeout in secs for initial handshake part
    timestampMaxDiff = 2  # maximum allowed difference in exchanged timestamps, in secs
    commonStartDelay = 3  # delay for determining the common start time

    # first stage: send first message ('helloszia')
    # and listen to the same message coming in from target address ("addr")
    startTime = time.time()
    while not successFlag and ((time.time() - startTime) < handshakeTimeout):
        ready_to_read, ready_to_write, in_error = select.select([socketComm],
                                                                [socketComm],
                                                                [socketComm],
                                                                0)
        # check if socket is in error state, return with error if yes
        if in_error:
            errorFlag = True
            print('Comm socket threw an error for select.select() during first stage of handshake!')
            return errorFlag, commonStartTime
        # if socket is empty, send init message
        elif not ready_to_read:
            if ready_to_write:
                socketComm.sendto('helloszia'.encode(), addr)
                print('Sending init message "helloszia" to remote comm socket at ' + str(addr))
        # if there is incoming packet, and it matches expected message + from right address, move on
        elif ready_to_read:
            incomingMess, incomingAddr = socketComm.recvfrom(512)
            print('Received message "' + incomingMess.decode() + '" from ' + str(incomingAddr))
            if bool(incomingMess == 'helloszia'.encode()) and incomingAddr == addr:
                print('Message matching init message, moving on...')
                socketComm.sendto('helloszia'.encode(), addr)
                successFlag = True
        # sleep between socket polls
        time.sleep(0.050)

    # check if timeout happened, return with error if yes
    if not successFlag:
        print('Handshake timed out during first stage!')
        errorFlag = True
        return errorFlag, commonStartTime

    # second stage: send timestamp and listen for incoming timestamp
    successFlag = False
    startTime = time.time()
    # for UDP packet we need a bytes array formed from the timestamp float
    timestampBytes = np.array(startTime).tobytes()
    while not successFlag and ((time.time() - startTime) < handshakeTimeout):
        ready_to_read, ready_to_write, in_error = select.select([socketComm],
                                                                [socketComm],
                                                                [socketComm],
                                                                0)
        # check if socket is in error state, return with error if yes
        if in_error:
            errorFlag = True
            print('Comm socket threw an error for select.select() during second stage of handshake!')
            return errorFlag, commonStartTime
        # if socket is empty, send timestamp in packet
        elif not ready_to_read:
            if ready_to_write:
                socketComm.sendto(timestampBytes, addr)
                print('Sending timestamp to remote comm socket at ' + str(addr))
        # if there is a packet in socket, check if it is a timestamp (float almost matching local timestamp)
        elif ready_to_read:
            incomingMess, incomingAddr = socketComm.recvfrom(512)
            if incomingMess != 'helloszia'.encode():
                incomingFloat = np.frombuffer(incomingMess)[0]
                print('Received message "' + str(incomingFloat) + '" from ' + str(incomingAddr))
                if abs(incomingFloat-startTime) > timestampMaxDiff:
                    print('Difference between timestamps in second stage of handshake is too large!')
                    return
                elif incomingAddr == addr:
                    print('Received valid timestamp from expected address, ready to start!')
                    socketComm.sendto(timestampBytes, addr)
                    commonStartTime = (startTime + incomingFloat)/2 + commonStartDelay
                    successFlag = True
        # sleep 10 ms between socket polls
        time.sleep(0.050)

    # check if timeout happened, return with error if yes
    if not successFlag:
        print('Handshake timed out during first stage!')
        errorFlag = True
        return errorFlag, commonStartTime
    # else simply return
    else:
        print('Handshake successful, negotiated shared start time!')
        return errorFlag, commonStartTime


def udpsender(socketOut, outAddr, startTime, packetNo=100, waitTime=0.01):
    """
    Function sending timestamps in UDP packets

    Mandatory inputs:
    socketOut :     Socket object, non--blocking UDP. Used for sending packets only.
    outAddr :       Tuple containing a network address (IPv4, port)
    startTime :     Float, timestamp in Unix time. The function waits with the first
                    socket state poll until "startTime"

    Optional inputs
    packetNo :      Numeric value, number of packets to send and expect. Defaults to 100.
    waitTime :      Numeric value, time between subsequent socket state polls, in secs.
                    Defaults to 0.01 (10 ms)

    No return vars, saves out timestamps ("outData").
    """

    # preallocate numpy arrays holding timestamps
    outData = np.zeros((packetNo, 1))

    # constants, counters
    roughWaitLimit = 0.1  # finish imprecise wait "roughWaitLimit" secs before target time

    # rough / imprecise wait till ("startTime"-"roughWaitLimit") with time.sleep()
    if time.time() >= startTime+roughWaitLimit:
        print('Shared start time is too close or already in the past, returning!')
        return
    else:
        time.sleep((startTime-time.time())-roughWaitLimit)

    # user message
    print('Starting UDP timestamp stream')

    # for loop for sending packets
    for packetOutIdx in range(packetNo):

        # for the first packet, wait precisely for "startTime"
        if packetOutIdx == 0:
            while time.time() < startTime:
                pass

        currentLoopStart = time.time()

        # poll socket states in a non-blocking way
        ready_to_read, ready_to_write, in_error = select.select([socketOut],
                                                                [socketOut],
                                                                [socketOut],
                                                                0)

        # if socket is in error, print message and return early
        if in_error:
            print('Socket for outgoing timestamps threw an error for select.select() at packetOutIdx ' + str(packetOutIdx) + ' !')
            return

        # write packet if we can, store sent timestamp in np.array
        if ready_to_write:
            outgoingTime = time.time()
            outgoingBytes = np.array(outgoingTime).tobytes()  # create bytes array from numpy float
            socketOut.sendto(outgoingBytes, outAddr)
            outData[packetOutIdx, 0] = outgoingTime

        # user message
        if packetOutIdx+1 % 50 == 0:
            print('Sent ' + str(packetOutIdx+1) + ' packets so far...')

        # wait before next iteration in while loop
        while time.time() < currentLoopStart+waitTime:
            pass

    # user end message
    print('Finished sending packets as planned. Packet counter at the end (started at zero):')
    print('Outgoing: ' + str(packetOutIdx))

    # save out timestamps
    np.savetxt('outData.csv', outData, delimiter=',')

    return


def udpreceiver(socketIn, packetNo=100):
    """
    Function receiving timestamps in UDP packets.

    Mandatory inputs
    socketIn :      Socket object, non-blocking UDP. Used for receiving packets only.

    Optional inputs
    packetNo :      Numeric value, number of packets to send and expect. Defaults to 100.

    Returns
    inData :        Numpy array, shaped [packetNo, 2], containing the two timestamps for each incoming packet:
                    (1) the one measured locally when packet was read out;
                    (2) the one in the packet, measured at send time at remote PC

    Saves out timestamps ("inData").
    """

    # preallocate numpy array holding timestamps
    inData = np.zeros((packetNo, 2))

    # packet counter
    packetInIdx = 0

    # user message
    print('Starting listening to incoming stream')

    while packetInIdx < packetNo:

        # poll socket states in a non-blocking way
        ready_to_read, ready_to_write, in_error = select.select([socketIn],
                                                                [socketIn],
                                                                [socketIn],
                                                                0)
        # if socket is in error, print message and return early
        if in_error:
            print('Socket for receiving timestamps threw an error for select.select() at packetInIdx ' + str(packetInIdx) + ' !')
            return

        # if there is anything to read, read packet, parse, store in np.array
        if ready_to_read:
            incomingTime = time.time()
            incomingBytes = socketIn.recv(128)
            incomingFloat = np.frombuffer(incomingBytes)[0]  # parsing the incoming packet as if it was a timestamp
            inData[packetInIdx, 0:2] = incomingTime, incomingFloat
            packetInIdx = packetInIdx + 1

        # user messages
        if packetInIdx+1 % 50 == 0:
            print('Received ' + str(packetInIdx+1) + ' packets so far...')

    # user end message
    print('Finished receiving packets as planned. Packet counter at the end (started at zero):')
    print('Incoming: ' + str(packetInIdx))

    # save out timestamps
    np.savetxt('inData.csv', inData, delimiter=',')

    return inData


def main(ip, portIn, portOut, portComm):
    '''
    Putting together the steps:
    (1) Create UDP sockets
    (2) Handshake - negotiate shared start time
    (3) Start UDP stream of timestamps + listen to incoming timestamps
    (4) Save data, report basic results
    (5) Close everything gracefully, return

    Inputs:
    ip :            Valid IPv4 address
    portIn :        Port number for incoming packets
    portOut :       Port number for outgoing packets
    portComm :      Separate port number for handshake procedure
    '''

    # create ports
    socketFlag, socketOut = openSocket(portOut)
    if not socketFlag:
        print('\n\nCould not create or bind UDP socket for portOut. Wtf.')
        sys.exit()
    socketFlag, socketIn = openSocket(portIn)
    if not socketFlag:
        print('\n\nCould not create or bind UDP socket for portIn. Wtf.')
        sys.exit()
    socketFlag, socketComm = openSocket(portComm)
    if not socketFlag:
        print('\n\nCould not create or bind UDP socket for portComm. Wtf.')
        sys.exit()

    # handshake
    print(socketComm)
    print(type(socketComm))
    print((ip, portComm))
    errorFlag, commonStartTime = handshake(socketComm, (ip, portComm))
    # check for errors
    if errorFlag:
        print('Ran into trouble while calling function "handshake", it returned an error.')
        print('Returned var commonStartTime is :' + str(commonStartTime))
        closeSockets([socketIn, socketOut, socketComm])
        return

    # start up a process sending packets
    packetSender = multiprocessing.Process(name='packetSender',
                                           target=udpsender,
                                           args=(socketOut,
                                                 (ip, portIn),
                                                 commonStartTime,))
    packetSender.start()

    # listen to incoming packets
    inData = udpreceiver(socketIn)

    # wait for packetSender to finish
    packetSender.join()

    # report basics
    diffsIn = inData[:, 0]-inData[:, 1]
    print('Median difference between remote and local timestamps for incoming packets: \n',
          "{:.8f}".format(np.median(diffsIn)))

    # user end message
    print('\n Its been fun, so long, etc. You take care now')
    closeSockets([socketIn, socketOut, socketComm])

    return


if __name__ == "__main__":
    # input arguments
    parser = argparse.ArgumentParser()
    # add arguments
    parser.add_argument(
        '-i',
        '--IP',
        nargs='?',
        type=str,
        default='localhost',
        help='Provide IP (ipv4) of remote PC in string.' +
             'Default = "localhost"')
    parser.add_argument(
        '-pi',
        '--portIn',
        nargs='?',
        type=int,
        default=9997,
        help='Port for incoming UDP packets.' +
             'Default = 9997')
    parser.add_argument(
        '-po',
        '--portOut',
        nargs='?',
        type=int,
        default=9998,
        help='Port for outgoing UDP packets.' +
             'Default = 9998')
    parser.add_argument(
        '-pc',
        '--portComm',
        nargs='?',
        type=int,
        default=9999,
        help='Port for handshake.' +
             'Default = 9999')
    # parse arguments
    args = parser.parse_args()

    main(args.IP,
         args.portIn,
         args.portOut,
         args.portComm)
