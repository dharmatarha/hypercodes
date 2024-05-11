"""

General version of behav_audio_merge.py.

USAGE: python3 audio_merge_gen.py DATA_DIR JSON_PATH

where
DATA_DIR:  Path to folder holding the audio files.
JSON_PATH: Path to fOut.json file holding the log file content for the audio files.

Saves out synched, stereo wavs into DATA_DIR, and also writes out the dataframe that holds all relevant data for each
audio file to "DATA_DIR/audio_merge_dataframe.pkl".


Tests:
>> get_tags_from_path('sub-S25_ses-dyad16_task-independent_desc-ReceivedAudio.wav')
{'ses': 'dyad16', 'task': 'independent', 'sub': 'S25', 'desc': 'ReceivedAudio'}
>> get_tags_from_path('some_folder1/some_folder2/sub-S25_ses-dyad16_task-independent_desc-ReceivedAudio.wav')
{'ses': 'dyad16', 'task': 'independent', 'sub': 'S25', 'desc': 'ReceivedAudio'}

"""


import json
import argparse
import wave
import sys
import numpy as np
from scipy.io import wavfile
from glob import glob
import os
import warnings
import pandas as pd


AUDIO_SAMPLING_RATE = 16000                              # Sampling rate of Joint Storytelling task audio files.
SEARCH_TERM = 'RecordedAudio.wav'                        # fOut.json content corresponds to RecordedAudio log files.
BYTES_PER_SAMPLE = 2                                     # 16 bit PCM encoding, thus 2 bytes = 1 sample
START_TIME_COMMON_POS = 0                                # The byte position in fOut.jon for each audio corresponding to startTimeCommon in the task
SECOND_TURN_START_TIME_POS = 1                           # The byte position in fOut.jon for each audio corresponding to the start of the second speech turn (and the end of the first turn)
JOINT_TASK_END_POS = 30                                  # The byte position in fOut.json for each audio corresponding to the end of the joint portion of the task
SPEECH_TURN_LENGTH = 30                                  # Speech turn length in seconds
DATAFRAME_PICKLE_FILENAME = 'audio_merge_dataframe.pkl'  # Filename for pickling the final dataframe into a file
STEREO_WAV_ENDING = 'synched_stereo'                     # Filename ending for stereo, synched wav files
SAMPLING_RATE_THRESH = 2.5                                 # Threshold for deviating from the nominal sampling rate in Hz, used when estimateing sampling rates from byte positions.


def get_tags_from_path(file_path):
    """
    Helper function for the Joint Storytelling task. Parses the file name for subject ID, dyad ID, task ID,
    and file descriptor if the file follows the naming convention:
    sub-[SUB_ID]_ses-[DYAD_ID]_task-[TASK_NAME]_desc-[FILE_DESCRIPTOR].*

    :param file_path:   String, path to a file or just the filename.
    :return: file_dict: Dictionary with key: value pairs: 'ses': [DYAD_ID], 'task': [TASK_NAME], 'sub': [SUB_ID],
                                                          'desc': [FILE_DESCRIPTOR]
    """

    # Get basename to parse.
    file_parts = os.path.basename(file_path).split('_')
    # Drop any file name extension from last part, simply by splitting at dot.
    file_parts[-1] = file_parts[-1].split('.')[0]

    # Get session (dyad) ID.
    ses_part = [p for p in file_parts if p.startswith('ses')]
    ses = ses_part[0].split('-')[1]
    # Get task ID.
    task_part = [p for p in file_parts if p.startswith('task')]
    task = task_part[0].split('-')[1]
    # Get subject ID.
    sub_part = [p for p in file_parts if p.startswith('sub')]
    sub = sub_part[0].split('-')[1]
    # Get file descriptor.
    desc_part = [p for p in file_parts if p.startswith('desc')]
    desc = desc_part[0].split('-')[1]

    # Check if the file contained the necessary session ID, task name, subject ID, and descriptor. Raise warning
    # if there is stg missing.
    if not ses or not task or not sub or not desc:
        warnings.warn('\nFile at ' + f + ' does not have session (dyad) ID or task name or subject ID or file descriptor!')

    file_dict = {'ses': ses,
                 'task': task,
                 'sub': sub,
                 'desc': desc}

    return file_dict


def audio_files_to_df(data_dir, pattern_in_filenames='RecordedAudio'):
    """
    Finds all files in "data_dir" that contain the pattern "pattern_in_filenames" in their names, and also have
    "sub", "ses", and "task" tags. The file paths are stored in a pandas dataframe where the values of the tags
    "sub", "ses", "task", and "desc" are also stored in separate columns.

    That is, for example, with an audio file "sub-99_ses-dyad99_task-control_RecordedAudio.wav" in "DATA_DIR",
    the command audio_files_to_df(DATA_DIR, pattern_in_filenames='RecordedAudio.wav')
    returns a 1-row pandas dataframe with columns "filepath", "filename", "sub", "ses", "task", and "desc" with values
    "DATA_DIR/sub-99_ses-dyad99_task-control_RecordedAudio.wav", "sub-99_ses-dyad99_task-control_RecordedAudio.wav",
    "99", "dyad99", "control", and "RecordedAudio", respectively.

    :param data_dir:              Path to dir containing the files we are searching for. Dir is globbed recursively.
    :param pattern_in_filenames:  String, some part of the filenames used for glob-bing recursively in "data_dir".
    :return: df:                  Pandas dataframe with one row for each file. Columns are "filepath", "filename",
                                  "session", "task", "subject", and "descriptor".
    """

    files_w_pattern = glob(os.path.join(data_dir, '**/*' + pattern_in_filenames + '*'), recursive=True)
    if not files_w_pattern:
        raise FileNotFoundError('Could not find any files matching the pattern!')

    # Init lists that will be the basis of the dataframe
    filepath_list = []
    filename_list = []
    ses_list = []
    task_list = []
    sub_list = []
    desc_list = []

    # Try to parse the filename for sub, ses / dyad, task and descriptor codes.
    for f in files_w_pattern:

        # Get session and subject ID, and task name from filename.
        tags_dict = get_tags_from_path(f)

        # Check if the filename contains session ID, task name, subject ID or descriptor. Raise warning if not.
        if not tags_dict['ses'] or not tags_dict['task'] or not tags_dict['sub'] or not tags_dict['desc']:
            warnings.warn('\nFile at ' + f + ' does not have session ID, task name, subject ID or descriptor!')

        # Else append values of tags to lists
        else:
            filepath_list.append(f)
            filename_list.append(os.path.basename(f))
            ses_list.append(tags_dict['ses'])
            task_list.append(tags_dict['task'])
            sub_list.append(tags_dict['sub'])
            desc_list.append(tags_dict['desc'])

    # Generate dataframe from lists
    df = pd.DataFrame({'filepath': filepath_list,
                       'filename': filename_list,
                       'session': ses_list,
                       'task': task_list,
                       'subject': sub_list,
                       'descriptor': desc_list})

    return df


def file_dict_to_df(df, file_info_dict, column_to_match='filename', new_column_name='byte_positions'):
    """"
    Helper function appending a pd dataframe with values in the dictionary "file_info_dict". The keys of
    "file_info_dict" must match the beginnings of strings in the "column_to_match" column of the dataframe to be paired
    uip. If there is a matching key, the corresponding value from the dictionary is added to the dataframe, in the
    variable "new_column_name".
    """
    df_new = df.copy()
    # Add variable "column_name" to dataframe if it is missing.
    if new_column_name not in df_new.columns:
        df_new.loc[:, new_column_name] = None
    # Loop through the dictionary keys, find matches in dataframe.
    for idx, key in enumerate(file_info_dict.keys()):
        index_list = [df_idx for df_idx, column_value in enumerate(df_new[column_to_match]) if column_value.startswith(key)]  # List of row indices.
        # If there was a match, check if it was unique, then store the dictionary value in the dataframe.
        if index_list:
            if len(index_list) == 1:
                print('Found match in data frame for dictionary key ' + key)
                print('Dataframe row idx: ' + str(index_list[0]) + ', matching value at ' +
                      column_to_match + ': ' + df_new.loc[index_list[0], column_to_match])
                df_new.at[index_list[0], new_column_name] = file_info_dict[key]
            else:
                warnings.warn('Multiple matches in dataframe for key ' + key)

    return df_new


def wavs_from_df(df, start_sample_no=START_TIME_COMMON_POS, end_sample_no=JOINT_TASK_END_POS,
                 wav_file_tag='synched_stereo', start_offset_s=0, end_offset_s=0, sr=AUDIO_SAMPLING_RATE):
    """
    Function to generate stereo audio wavs from the mono wav files, trimmed to the start and end of the joint task.
    Input is a pandas dataframe, where each row corresponds to the data of one mono audio file, with variables on
    the session, subject, task, fOut logfile content, etc. The outputs are the wav files themselves, with their paths
    added to the dataframe in column "synched_wav_path".

    Input dataframe must have columns

    :param df:              Pandas dataframe.
    :param start_sample_no: Which audio sample position in df column "audio_samples" corresponds to joint task start.
    :param end_sample_no:   Which audio sample position in df column "audio_samples" corresponds to joint task end.
    :param wav_file_tag:    String, added to the base filenames of the original audio files to generate output
                            wav names.
    :param start_offset_s:  Numeric value, offset from the sample no defined in arg start_sample_no, in seconds. Offset
                            is applied at trimming. Defaults to 0.
    :param end_offset_s:    Numeric value, offset from the sample no defined in arg end_sample_no, in seconds. Offset
                            is applied at trimming. Defaults to 0.
    :param sr:              Numeric value, sampling rate in HZ, required for calculating offset if args start_offset_s or
                            end_offset_s are supplied. Defaults to module constant AUDIO_SAMPLING_RATE.
    :return: df_new:        Pandas dataframe, same as "df", but with extra column ("synched_wav_path") holding the path
                            of the output wav files.
    """
    df_new = df.copy()
    df_new.loc[:, 'synched_wav_path'] = None

    # Loop through audio files
    for row_idx in df_new.index:

        # Only take the row into account if the corresponding stereo has not been written out yet.
        if not df_new.loc[row_idx, 'synched_wav_path']:

            # Is the corresponding other mono wav file available (in the dataframe)?
            mask = (df_new.loc[:, 'session'] == df_new.loc[row_idx, 'session']) & \
                   (df_new.loc[:, 'task'] == df_new.loc[row_idx, 'task']) & \
                   (df_new.loc[:, 'descriptor'] == df_new.loc[row_idx, 'descriptor']) & \
                   (df_new.loc[:, 'subject'] != df_new.loc[row_idx, 'subject'])
            # Mask should be a boolean pd series with exactly one match.
            if sum(mask) != 1:
                warnings.warn('Could not find matching pair for dataframe row ' + str(row_idx + '!!!'))

            # Select both wav files and the start and end samples.
            # First subject:
            wav1_path = df_new.loc[row_idx, 'filepath']
            samples1 = df_new.loc[row_idx, 'audio_samples'][[start_sample_no, end_sample_no]]
            # Second subject
            wav2_path = df_new.loc[mask, 'filepath'].values[0]  # Workaround, as the masking returns a pd Series object.
            samples2 = df_new.loc[mask, 'audio_samples'].values[0]  # Workaround, as the masking returns a pd Series object.
            samples2 = samples2[[start_sample_no, end_sample_no]]

            # Apply offsets if needed
            if start_offset_s:
                samples1[0] = samples1[0] + start_offset_s * sr
                samples2[0] = samples2[0] + start_offset_s * sr
            if end_offset_s:
                samples1[1] = samples1[1] + end_offset_s * sr
                samples2[1] = samples2[1] + end_offset_s * sr

            # User feedback
            print('\nFound valid wav files for session: ' + df_new.loc[row_idx, 'session'] + ', task: ' +
                  df_new.loc[row_idx, 'task'] + ':')
            print(wav1_path)
            print(wav2_path)

            # Load wavs, trim them
            # First subject:
            frate1, channel1 = wavfile.read(wav1_path)
            channel1 = channel1[samples1[0]:samples1[1]]
            # Second subject:
            frate2, channel2 = wavfile.read(wav2_path)
            channel2 = channel2[samples2[0]:samples2[1]]

            # Sanity check
            if frate1 != sr or frate2 != sr:
                raise ValueError('Wrong sampling rates!!!')

            # Further trim to the shorter wav, there is usually some small difference.
            minsize = min(channel1.size, channel2.size)
            channel1 = channel1[:minsize]
            channel2 = channel2[:minsize]
            # Create matrix for stereo output wav.
            channels = np.vstack([channel1, channel2])

            # Write stereo output wav
            stereo_filename = '_'.join(['ses-' + df_new.loc[row_idx, 'session'],
                                        'task-' + df_new.loc[row_idx, 'task'],
                                        'desc-' + df_new.loc[row_idx, 'descriptor'],
                                        wav_file_tag + '.wav'])
            stereo_dir = os.path.split(wav1_path)[0]
            stereo_path = os.path.join(stereo_dir, stereo_filename)
            print('Writing stereo wav file to ' + stereo_path)
            print(channels.shape)
            wavfile.write(stereo_path, AUDIO_SAMPLING_RATE, channels.T)

            # Update dataframe
            df_new.loc[row_idx, 'synched_wav_path'] = stereo_path
            df_new.loc[mask, 'synched_wav_path'] = stereo_path

    return df_new


def sampling_rates_in_df(df, nominal_sr=AUDIO_SAMPLING_RATE,
                         byte_position_indices=(SECOND_TURN_START_TIME_POS, JOINT_TASK_END_POS),
                         turn_length_s=SPEECH_TURN_LENGTH,
                         sampling_rate_thresh=SAMPLING_RATE_THRESH):
    """
    Function to loop over all wav data in dataframe "df" to check if their estimated sampling rates matches the nominal
    sampling rate. We assume that elapsed time between byte_position_indices is
    (byte_position_indices[1] - byte_position_indices[0]) * turn_length_s, that is, byte_positions correspond to speech
    turn starts (and endings) and not other events.

    :param df:                    Pandas dataframe.
    :param nominal_sr:            Expected / nominal sampling rate in Hz.
    :param byte_position_indices: Tuple of indices to select values from the "byte_positions" df column with known
                                  timing to pass to estimate_sampling_rate for each wav.
    :param turn_length_s:         Speech turn length in seconds, we assume that elapsed time between byte_position_indices is
    :param sampling_rate_thresh:  Threshold for providing a warning about possibly problematic sampling rates.
    :return: df_new:              Pandas dataframe, same as "df", but with extra column ("estimated_sr") holding the
                                  estimated sampling rate value (in Hz).
    """
    df_new = df.copy()
    df_new.loc[:, 'estimated_sr'] = None

    # Nominal elapsed time between the byte positions of the audio as indexed by byte_position_indices.
    elapsed_time = (byte_position_indices[1] - byte_position_indices[0]) * turn_length_s

    # Loop through audio files
    for row_idx in df_new.index:

        byte_pos = df_new.loc[row_idx, 'byte_positions']
        byte_position_tuple = (byte_pos[byte_position_indices[0]], byte_pos[byte_position_indices[1]])

        sr_estimate = estimate_sampling_rate(byte_position_tuple, elapsed_time)

        df_new.loc[row_idx, 'estimated_sr'] = sr_estimate

        if abs(nominal_sr - sr_estimate) > sampling_rate_thresh:
            warnings.warn('\nPotentially bad sampling rate (' + str(sr_estimate) + ') at row ' + str(row_idx) + '!')
            print('Filename:', df_new.loc[row_idx, 'filename'])

    return df_new


def estimate_sampling_rate(byte_positions, elapsed_time_s, bytes_per_sample=BYTES_PER_SAMPLE):
    """
    Helper to estimate the sampling rate from known byte positions and timing information from an audio recording.
    :param byte_positions:   Tuple with two elements, starting and ending byte positions for a segment with known timing.
    :param nominal_time_s:   Elapsed time between byte_positions[0] and byte_positions[1] in seconds.
    :param bytes_per_sample: Number of bytes per audio sample, defaults to module constant BYTES_PER_SAMPLE.
    :return:                 Estimated sampling rate in Hz.
    """
    return (byte_positions[1] - byte_positions[0]) / (elapsed_time_s * bytes_per_sample)



def main():
    # Parse input arguments.
    parser = argparse.ArgumentParser()
    # Input arg "info_json" points to fOut.json file holding timing info for all audio files.
    parser.add_argument('data_dir', type=str, help='Path to folder holding audio (wav) files.')
    parser.add_argument('json_path', type=str, help='Path to fOut.json file.')
    # Parse arguments, get list.
    args = parser.parse_args()

    print('\nCalled audio_merge_gen with args:')
    print('Data dir:', args.data_dir)
    print('JSON path:', args.json_path)

    # Read in json with timestamps.
    with open(args.json_path, 'r') as json_file:
        json_data = json.load(json_file)
    print('\nLoaded json into dictionary with length ' + str(len(json_data)) + '.')

    # List all files matching SEARCH_TERM and store their details in a dataframe.
    df = audio_files_to_df(args.data_dir, pattern_in_filenames=SEARCH_TERM)
    print('\nListed audio files from data dir and listed them in a dataframe (' + str(len(df)) + ' files total).')
    print('Dataframe variables are ', list(df.columns))

    # Match timestamps from json to audio files.
    print('\nMatching json content to audio file, appending dataframe rows...')
    df = file_dict_to_df(df, json_data, column_to_match='filename', new_column_name='byte_positions')
    print('Done. Dataframe variables are ', list(df.columns))

    # Turn byte positions into samples and time
    print('\nTransforming byte positions into audio samples and timestamps...')
    df.loc[:, 'audio_samples'] = None
    df.loc[:, 'timestamps'] = None
    for row_idx in df.index:
        # Turn byte_positions list into numpy array
        df.at[row_idx, 'byte_positions'] = np.asarray(df.loc[row_idx, 'byte_positions'])
        # Get samples from byte_positions, as array
        df.at[row_idx, 'audio_samples'] = (df.loc[row_idx, 'byte_positions'] / BYTES_PER_SAMPLE).astype(int)
        # Get timestamps from samples
        df.at[row_idx, 'timestamps'] = df.loc[row_idx, 'audio_samples'] / AUDIO_SAMPLING_RATE
    print('Done. Dataframe variables are ', list(df.columns))

    # Check for bad sampling rates
    print('\nChecking for bad sampling rates...')
    df = sampling_rates_in_df(df, byte_position_indices=(SECOND_TURN_START_TIME_POS, JOINT_TASK_END_POS))
    print('Done. Dataframe variables are ', list(df.columns))

    # Define start and stop points for JointStory portion of the audio, edit wav files accordingly
    print('\nPairing up wav files, trimming, and writing them out in stereo...')
    df = wavs_from_df(df, start_sample_no=SECOND_TURN_START_TIME_POS,
                      end_sample_no=JOINT_TASK_END_POS,
                      wav_file_tag=STEREO_WAV_ENDING,
                      start_offset_s=-1 * SPEECH_TURN_LENGTH)
    print('Done. Dataframe variables are ', list(df.columns))

    # Save out dataframe for troubleshooting, use pickle
    df_filename = os.path.join(args.data_dir, DATAFRAME_PICKLE_FILENAME)
    df.to_pickle(df_filename)
    print('\nSaved out the dataframe out with pickle to ' + df_filename)


if __name__ == '__main__':
    main()
