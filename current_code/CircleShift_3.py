# necessary imports
import numpy as np
import random


def Circle_Shift_3(data_matrix, shift_point=float("nan"), shift_axis=1, rand_shift=True, shuffleRange=0):
    """
    Circle shifts an input data matrix. Essentially takes the bottom n rows of a matrix and tosses them at the top.
    The variable n is determined by the shift_point variable (or, to randomly select it with each call,
    set rand_shift to true).

    This approach is a) suggested for use in ISC analyses by Nastase et al., 2019,
    b) used in an Intersubject Correlations tutorial by Juha Lahnakoski and Luke Chang, and
    c) described in section 7.4 of Lancaster et al., 2008.

    Input requirements:
    -Matrix: 2D np array with rows = samples and columns = variables.
    -shift_point: int between 0 and the max number of samples (only matters if rand_shift = False)
    -shift_axis: either 0 or 1 (ints)
    -rand_shift: Boolean (True or False)
    -shuffleRange: int indicating number of TRs over which to implement the circle shifting process.
               The default of 0 ignores the segmentation process and implements circle shifting over the
               entire time series.
    -numpy loaded as np
    -random loaded as random

    Input suggestions:
    -If rand_shift is set to true, shift points will be automatically picked at random with each function call.
     Points are picked within the range of (0,number of samples) with the endpoints excluded
    -If rand_shift is set to false, you must supply the shift point via the shift point variable.
    -The matrix is automatically shifted along columns (i.e. the bottom n rows are moved to the top).
     If your data has rows = variables and columns = samples, then set shift_axis = 0.

    Returns:
    -A shifted version of data_matrix with the bottom n samples put at the start
    """

    ### Preprocessing and input checks ###

    # check those inputs
    if not isinstance(data_matrix, np.ndarray) or len(data_matrix.shape) != 2:
        print('Input arg "data_matrix" should be a 2D numpy array of reals!')
        return 0
    if shift_axis not in [0, 1]:
        print('Input arg "shift_axis" should be 0 or 1!')
        return 0
    if rand_shift not in [True, False]:
        print('Input arg "rand_shift" should be True or False!')
        return 0
    if not isinstance(shuffleRange, int) or shuffleRange < 0 or shuffleRange == 1:
        print('Input argument "shuffleRange" must be a positive integer greater than or equal to 2!')
        return 0

    # transpose data_matrix if shift_axis = 0 -> rows will now be samples, cols = vars
    if shift_axis == 0:
        data_matrix = np.transpose(data_matrix)

    # if shuffleRange is zero, set the shuffling range to the total number of samples in the time series
    if shuffleRange == 0:
        shuffleRange = data_matrix.shape[0]

    # if shuffleRange does not go into the number of TRs evenly...
    if not (data_matrix.shape[0] / shuffleRange) % 1 == 0:

        # feedback
        print('WARNING! The shuffling range does not go into the number of TRs evenly!')

        # figure out an uneven segmentation scheme
        numMainSegs = np.floor(data_matrix.shape[0] / shuffleRange)
        segments = numMainSegs + 1
        mainSeg = shuffleRange
        lastSeg = data_matrix.shape[0] - mainSeg * numMainSegs

        # feedback
        print('Circle shifting will occur across ' + str(numMainSegs) + ' segments of ' + str(mainSeg) + ' samples and one segment of ' + str(lastSeg) + ' samples.')

        # get segment sample indices
        segs = [[]] * segments
        for SEG in range(segments):

            if SEG == 0:
                segs[SEG] = np.arange(mainSeg).astype(int)
            elif SEG == numMainSegs:
                segs[SEG] = (segs[SEG - 1][np.arange(lastSeg).astype(int)] + mainSeg).astype(int)
            else:
                segs[SEG] = (segs[SEG - 1] + mainSeg).astype(int)

    else:

        # get segment sample indices
        mainSeg = shuffleRange
        segments = int(data_matrix.shape[0] / shuffleRange)
        segs = [[]] * segments
        for SEG in range(segments):
            if SEG == 0:
                segs[SEG] = np.arange(mainSeg).astype(int)
            else:
                segs[SEG] = (segs[SEG - 1] + mainSeg).astype(int)

    # preallocate the shifted data matrix
    shifted_data_matrix = np.empty(data_matrix.shape)

    # for each segment...
    for SEG in range(len(segs)):

        # get current segment of inputted data matrix
        tmp = data_matrix[segs[SEG],:]

        # select a random shift value in the range (0,nsamples) if rand_shift = True.
        if rand_shift == True:
            shift_point = random.randint(1, np.size(tmp, 0) - 1)

        # If rand_shift = False, make sure shift_point is an int in the specified range
        elif isinstance(shift_point, int) == False or shift_point <= 0 or shift_point >= np.size(tmp, 0):
            print('Input arg "shift_point" must be an int in the range of (0,number of samples) - endpoints excluded')
            return 0

        # take the bottom n rows and move them to the top
        shifted_data_matrix[segs[SEG],:] = np.vstack((tmp[shift_point:], tmp[0:shift_point]))

    # return data in desired format if axis switched
    if shift_axis == 0:
        shifted_data_matrix = np.transpose(shifted_data_matrix)

    return shifted_data_matrix