#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Script for google's new speech-to-text api, using long_running_recognition.

Notes:
- Only one alternative is queried (max_alternatives of the config is not set)
- You need to edit the path for your own credentials file
- Only files in a bucket in google cloud storage can be used


    MANDATORY INPUT:

AUDIO: Path to an object (file) in google storage, including the bucket
       and object names. E.g.: audio = "my_bucket/my_audio_file".
       The full URI (request endpoint) will then look like this:
       "gs://my_bucket/my_audio_file",
       and that's enough for the speech-to-text app.
       Using a bucket allows for long audio (> 1 min).
       WARNING: You need to know in advance the encoding / sample rate / etc
       of the audio file, as here we only point to it. If you are not using
       a default LINEAR16 encoded, 16 kHz wav file, set the appropriate
       optional arguments.


    OPTIONAL INPUT:

CONTEXT: Speech context. Phrases that could help the recognition, especially
         helpful with idiosyncratic names. Google expects a list of phrases
         for the speech context object, where each phrase is max 100 chars,
         and there are no more than 500 phrases alltogether (there is also
         and overall cap of 10 000 chars though).
         CONTEXT is a text file, with each line as one phrase.
         Default is no context.

LANGUAGE: Speech-to-text can do many different languages now and also
          many dialects. Specifying the dialect can help the transcription.
          LANGUAGE should be a language code (string).
          For a list of supported language codes see:
          https://cloud.google.com/speech/docs/languages
          Default is 'en-US'

ENCODING: Audio file encoding, we need to supply this to the
          speech-to-text engine.
          String, defaults to 'LINEAR16'.
          For a list of supported values see:
          https://cloud.google.com/speech/reference/rest/v1beta1/RecognitionConfig#AudioEncoding

SAMPRATE: Audio file sampling rate, we need to supply this to the
          speech-to-text engine.
          Int, defaults to 16000.
          Values between 8000-48000 are supported, but 16000 is advised.


    OUTPUTS:

TRANSCRIPT FILE: txt file containing the transcript. No punctuation.

WORDLIST FILE: csv file listing all the words. Each row has the format:
               [startTime, endTime, word]
               WARNING: endTime is unreliable. It always equals the
               startTime of the next word, meaning that pauses count
               into the legth of the last word


Created on Sat Jan 20 15:29:19 2018

@author: adam.boncz@gmail.com
"""


#%% IMPORTS

import argparse
import csv
import time
import os
from google.cloud import speech


#%% Set GOOGLE_APPLICATION_CREDENTIALS environment variable

def setCredentials():
    credPath = '/home/adamb/Documents/Python/google_serv_Account/\
SpeechToText-e9f958389223.json'
    os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = credPath


#%% Write out results

# Response object levels seem to be:
#
# response,
#   result in response.results,
#     alternative in result.alternatives,
#     transcript in alternative.transcript,
#     confidence in alternative.confidence
#     word in alternative.words
#       start_time in word.start_time
#       end_time in word.end_time
#       word in word.word
#
# Write out one file with only the trancripts,
# and another one with timing data for each word


def writeResponse(saveFileBase, response):

    # First write out the pure transcript parts
    saveFile = saveFileBase + '_transcript.txt'
    with open(saveFile, 'w') as f:
        for part in range(len(response.results)):
            f.write(response.results[part].alternatives[0].transcript +
                    ' ')

    # write out each word with start and end times into a csv,
    # first column: start time, second: end time, third: word
    saveFile = saveFileBase + '_words.csv'
    with open(saveFile, 'w') as file:
        writer = csv.writer(file, delimiter=',')
        # iterate over all parts of results
        for part in range(len(response.results)):
            # iterate over all words
            for word in range(len(response.results[part].alternatives[0].words)):
                # derive start and end times as floats,
                # plus the word itself as string
                wordInfo = response.results[part].alternatives[0].words[word]
                startTime = float(str(wordInfo.start_time.seconds) +
                                  '.' + str(wordInfo.start_time.nanos))
                endTime = float(str(wordInfo.end_time.seconds) +
                                '.' + str(wordInfo.end_time.nanos))
                currentWord = wordInfo.word
                writer.writerow([startTime, endTime, currentWord])


#%% Check input arguments

def checkInputs(args):

    # Check speechcontext
    if args.speechcontext:
        # if file, try to read it in, and create a list of phrases
        if os.path.exists(args.speechcontext):
            with open(args.speechcontext, 'r') as speechFile:
                context = [line.strip() for line in speechFile.readlines()]
                print('\nLoaded speech context file, found',
                      len(context), 'phrases')
        else:
            print('\nFound no file with the name given for speechcontext')
    # if there was no speechcontext argument, or was empty
    else:
        print('\nReceived no speechcontext, '
              'we trust google to do magic on its own')
        context = []

    # Check language
    if not args.language:
        language = 'en-US'
    else:
        language = args.language
    print('\nLanguage code is set to ' + language)

    # Check encoding
    if not args.encoding:
        encoding = 'LINEAR16'
    else:
        encoding = args.encoding
    print('\nEncoding is set to ' + encoding)

    # Check sampling rate
    if not args.samprate:
        samprate = 16000
    else:
        samprate = args.samprate
    print('\nSampling rate is set to', samprate)

    return context, language, encoding, samprate


#%% Use speech-to-text service

def transcribe(audio, context, language, encoding, samprate):

    # set start time
    startTime = time.time()

    # init service
    client = speech.SpeechClient()
    # use long_running_recognize, aka asynchronous service,
    # supply parameters in audio and config objects
    operation = client.long_running_recognize(
                audio=speech.types.RecognitionAudio(
                      uri='gs://'+audio
                      ),
                config=speech.types.RecognitionConfig(
                       encoding=encoding,
                       language_code=language,
                       sample_rate_hertz=samprate,
                       enable_word_time_offsets=True,
                       speech_contexts=[speech.types.SpeechContext(
                                       phrases=context)]
                       ),
                )

    # results
    response = operation.result()
#    print(results)

    # feedback
    print('\nReceived response, transcription took '
          '{0:.2f} secs'.format(time.time()-startTime))

    return response


#%% Main

def main():

    # Argument parsing
    parser = argparse.ArgumentParser(description='Script for connecting '
                                     'to Google\'s speech-to-text service. '
                                     'Read the docstring in the script.')
    # Adding arguments. 'Audio' is mandatory, 'speechcontext', 'language',
    # 'encoding' and 'samprate' are optional
    # Audio
    parser.add_argument(
        'audio',
        type=str,
        help='URI to audio file in google storage '
             '(e.g. \'my_bucket/my_audio_object\'.')
    # Speechcontext
    parser.add_argument(
        '-s',
        '--speechcontext',
        type=str,
        const=None,
        default=None,
        help='Text file with the speechcontext,'
             'each line is a phrase supplied to the engine. '
             'Defaults to empty list')
    # Languagehttps://hynek.me/articles/serialization/
    parser.add_argument(
        '-l',
        '--language',
        type=str,
        const=None,
        default=None,
        help='Language/dialect code supplied to the engine. '
             'Defaults to \'en-US\'')
    # Encoding
    parser.add_argument(
        '-e',
        '--encoding',
        type=str,
        const=None,
        default=None,
        help='Audo encoding type supplied to the engine. '
             'Defaults to \'LINEAR16\'')
    # Sampling rate
    parser.add_argument(
        '-r',
        '--samprate',
        type=int,
        const=None,
        default=None,
        help='Audo sampling rate supplied to the engine. '
             'Defaults to 16000')

    # parse arguments
    args = parser.parse_args()

    # check inputs
    context, language, encoding, samprate = checkInputs(args)

    # set env variable for authentication
    setCredentials()

    # use speech-to-text
    response = transcribe(args.audio, context, language, encoding, samprate)

    # write out results,
    # first create saveFile name
    saveFileBase = os.path.splitext(os.path.split(args.audio)[1])[0]
    writeResponse(saveFileBase, response)


#%% GO

if __name__ == '__main__':
    main()
