import numpy as np
import matplotlib.pyplot as mplot
#from numpy.random import default_rng #***jd edit***
from numpy.random import random_sample #***jd edit***
from math import pi

def naiveColumnCorr(a, b):
    """
    Naive, slow, "baseline" function correlating corresponding columns of two matrices.
    Uses a for loop across columns.
    """
    c = np.zeros((a.shape[1]))
    for i in range(a.shape[1]):
        c[i] = np.corrcoef(a[:, i], b[:, i])[0, 1]

    return c


def fastColumnCorr(a, b):
    """
    Fast function for correlating corresponding columns of two matrices.
    Uses numpy.einsum to avoid loops and do computations directly on matrices.
    About ~ 10 times faster than the naive approach in 'naiveColumnCorr'.
    Inputs are 2D numpy arrays with the same shape, both sized samples X vars.
    NOTES:
    Could be further optimized using numpy.einsum_path for contraction order before first use,
    then simply calling einsum with that order subsequently. However, it only seems to give a
    few percents at best.
    contr_order = np.einsum_path("ij,ij->j", aa, bb, optimize='optimal')
    cov = np.einsum("ij,ij->j", aa, bb, optimize=contr_order[1])
    """
    # subtract the means from each var, in both matrices
    aa = a - (np.sum(a, 0) / a.shape[0]) # compute a - mean(a)
    bb = b - (np.sum(b, 0) / b.shape[0]) # compute b - mean(b)

    # multiply and sum across rows, that is, get dot products of column pairs
    cov = np.einsum("ij,ij->j", aa, bb)

    # for normalization we need the variances, separately for each var
    var_a = np.sum(aa ** 2, 0)
    var_b = np.sum(bb ** 2, 0)

    return cov / np.sqrt(var_a*var_b)


def phase_scrambling(data_matrix, fft_axis=0):
    """
    FOR REAL DATA ONLY, NOT COMPLEX!
    Phase-scrambling function for matrices. Preserves the original covariance structure.
    After FFT, we add a random phase vector to the FFT components of all time series / vars and do inverse FFT,
    as described in Prichard and Theiler (1994, Generating surrogate data for time series with several
    simultaneously measured variables. Physical review letters, 73(7), 951).
    The returned phase-scrambled data has the same power spectrum as the original but is linearly independent
    (zero expected correlation). Covariance structure is preserved, meaning that linear dependencies
    across time series is the same in the phase-scrambled data as in the original.
    Inputs:
    data_matrix:        2D numpy array of reals. Time series (Vars) X samples by default,
                        set fft_axis if samples X time series.
    fft_axis:           Axis along which FFT / iFFT is calculated. Defaults to 0,
                        meaning that FFT is calculated across rows (= each column is a separate time series / var)
    Outputs:
    data_scrambled:     2D numpy array of reals, contains the phase-scrambled data.
                        Same size and dimensions as input "data_matrix".
                        Returns 0 if input checks failed.
    TODO:
    - look into implementation with FFTW, which is supposedly faster with repetitive usage (our use case)
    """

    # input checks
    if not isinstance(data_matrix, np.ndarray) or len(data_matrix.shape) != 2 or np.iscomplex(data_matrix).any():
        print('Input arg "data_matrix" should be a 2D numpy array of reals!')
        return 0
    if fft_axis not in [0, 1]:
        print('Input arg "fft_axis" should be 0 or 1!')
        return 0

    # if fft_axis != 0, transpose the data
    transposeFlag = False
    if fft_axis == 1:
        data_matrix = np.transpose(data_matrix)
        transposeFlag = True

    # do forward FFT, use version for reals, treat data as if vars were in columns
    data_fft = np.fft.rfft(data_matrix, axis=0)

    # convert to polar coordinates (amplitude/magnitude + phase)
    data_fft_amp = np.abs(data_fft)
    data_fft_angle = np.angle(data_fft)

    # get random phase vector  (values between 0 - 2pi) for all FFT components that are not real by definition
    #rng = default_rng()  # new recommended method for random values
    if data_matrix.shape[0] % 2 == 0:  # if even, first and last components are real
        rand_phases = np.hstack(([0] ,(np.random.random_sample((data_fft.shape[0]-2)) * 2 * pi), [0])) #***jd edit***
    else:  # otherwise only the first component is real
        rand_phases = np.hstack(([0], (np.random.random_sample((data_fft.shape[0] - 1)) * 2 * pi))) #***jd edit***

    # add random phases to the angles of FFT components of all time series / vars,
    # addition is with broadcasting (newaxis is needed for broadcasting)
    data_fft_angle_rand = data_fft_angle + rand_phases[:, np.newaxis]

    # transform back from polar to cartesian, using the randomized phases but the original magnitude / amplitude values
    data_scrambled_fft = data_fft_amp * np.exp(1j * data_fft_angle_rand)  # returns complex FFT coefficients

    # do inverse FFT
    data_scrambled = np.fft.irfft(data_scrambled_fft, n=data_matrix.shape[0], axis=0)

    # transpose if necessary
    if transposeFlag:
        data_scrambled = np.transpose(data_scrambled)

    return data_scrambled

def phase_scrambling_2(data_matrix, fft_axis=0, shuffleRange=0):
    """
    FOR REAL DATA ONLY, NOT COMPLEX!
    Phase-scrambling function for matrices. Preserves the original covariance structure.
    After FFT, we add a random phase vector to the FFT components of all time series / vars and do inverse FFT,
    as described in Prichard and Theiler (1994, Generating surrogate data for time series with several
    simultaneously measured variables. Physical review letters, 73(7), 951).
    The returned phase-scrambled data has the same power spectrum as the original but is linearly independent
    (zero expected correlation). Covariance structure is preserved, meaning that linear dependencies
    across time series is the same in the phase-scrambled data as in the original.
    Inputs:
    data_matrix:            2D numpy array of reals. Time series (Vars) X samples by default,
                            set fft_axis if samples X time series.
    fft_axis:               Axis along which FFT / iFFT is calculated. Defaults to 0,
                            meaning that FFT is calculated across rows (= each column is a separate time series / var)
    shuffleRange:           Int indicating number of TRs over which to implement the phase scrambling process.
                            The default of 0 ignores the segmentation process and implements circle shifting over the
                            entire time series.
    Outputs:
    scrambled_data_matrix:  2D numpy array of reals, contains the phase-scrambled data.
                            Same size and dimensions as input "data_matrix".
                            Returns 0 if input checks failed.
    TODO:
    - look into implementation with FFTW, which is supposedly faster with repetitive usage (our use case)
    """

    # input checks
    if not isinstance(data_matrix, np.ndarray) or len(data_matrix.shape) != 2 or np.iscomplex(data_matrix).any():
        print('Input arg "data_matrix" should be a 2D numpy array of reals!')
        return 0
    if fft_axis not in [0, 1]:
        print('Input arg "fft_axis" should be 0 or 1!')
        return 0
    if not isinstance(shuffleRange, int) or shuffleRange < 0 or shuffleRange == 1:
        print('Input argument "shuffleRange" must be a positive integer greater than or equal to 2!')
        return 0

    # if fft_axis != 0, transpose the data
    transposeFlag = False
    if fft_axis == 1:
        data_matrix = np.transpose(data_matrix)
        transposeFlag = True

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
        print('Circle shifting will occur across ' + str(numMainSegs) + ' segments of ' + str(
            mainSeg) + ' samples and one segment of ' + str(lastSeg) + ' samples.')

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
    scrambled_data_matrix = np.empty(data_matrix.shape)

    # for each segment...
    for SEG in range(len(segs)):

        # do forward FFT, use version for reals, treat data as if vars were in columns
        data_fft = np.fft.rfft(data_matrix[segs[SEG],:], axis=0)

        # convert to polar coordinates (amplitude/magnitude + phase)
        data_fft_amp = np.abs(data_fft)
        data_fft_angle = np.angle(data_fft)

        # get random phase vector  (values between 0 - 2pi) for all FFT components that are not real by definition
        # rng = default_rng()  # new recommended method for random values
        if data_fft.shape[0] % 2 == 0:  # if even, first and last components are real
            rand_phases = np.hstack(([0] ,(np.random.random_sample((data_fft.shape[0]-2)) * 2 * pi), [0])) # ***jd edit***
        else:  # otherwise only the first component is real
            rand_phases = np.hstack(([0], (np.random.random_sample((data_fft.shape[0] - 1)) * 2 * pi))) # ***jd edit***

        # add random phases to the angles of FFT components of all time series / vars,
        # addition is with broadcasting (newaxis is needed for broadcasting)
        data_fft_angle_rand = data_fft_angle + rand_phases[:, np.newaxis]

        # transform back from polar to cartesian, using the randomized phases but the original magnitude / amplitude values
        data_scrambled_fft = data_fft_amp * np.exp(1j * data_fft_angle_rand)  # returns complex FFT coefficients

        # do inverse FFT
        data_scrambled = np.fft.irfft(data_scrambled_fft, n=data_matrix[segs[SEG],:].shape[0], axis=0)

        # transpose if necessary
        if transposeFlag:
            data_scrambled = np.transpose(data_scrambled)

        scrambled_data_matrix[segs[SEG],:] = data_scrambled

    return scrambled_data_matrix


def phase_scrambling_tests(data_matrix, data_scrambled, fft_axis=0, epsilon=1e-10):
    """
    Tests for the phase_scrambling function:
    (1) compare original FFT amplitudes to scrambled data FFT amplitudes
    (2) compare original covariance matrix to scrambled data covariance matrix
    (3) check if correlations between original and corresponding scrambled time series are around 0
    In the third test, we expect the correlation coefficients to show a normal
    distribution around 0, with a "small" std. To keep things simple, we do not fit
    a normal distribution or try a formal statistical test, but plot the histogram
    of the values and decide the test on the basis of the mean and median values
    (we check if they are "close" to zero, meaning < 0.05).
    IMPORTANT: For the second check we calculate the covariance matrices, so for really large data
    (e.g. tens of thousands of variables) consider the memory requirements of that step
    (~ 1.6 GB for 10^3 variables, considering we need two matrices). The function does not have
    internal checks for that.
    Inputs:
    data_matrix:        2D numpy array of reals. Original data set before phase scrambling.
                        Time series (Vars) X samples by default, set fft_axis if samples X time series.
    data_scrambled:     2D numpy array of reals, phase scrambled version of "data_matrix".
                        Same size and shape as "data_matrix".
    fft_axis:           Axis along which FFT / iFFT is calculated. Defaults to 0,
                        meaning that FFT is calculated across rows (= each column is a separate time series / var).
    epsilon:            Numeric value, threshold for machine accuracy. Tests 1 and 2 are considered "passed"
                        (that is, output "test_results" values set to True), if numeric inaccuracies
                        are below the threshold "epsilon". Defaults to 1e-10.
    Output:
    test_results:       List of booleans, 3-element long. Each boolean value corresponds
                        to pass (True) / fail (False) on a test.
                        The three values correspond to the (1) FFT amplitude test, (2) covariance matrix test,
                        and (3) correlations test.
                        Returns 0 if input checks failed.
    """

    # input checks
    if not isinstance(data_matrix, np.ndarray) or len(data_matrix.shape) != 2 or np.iscomplex(data_matrix).any():
        print('Input arg "data_matrix" should be a 2D numpy array of reals!')
        return 0
    if type(data_matrix) != type(data_scrambled) or data_matrix.shape != data_scrambled.shape or np.iscomplex(data_scrambled).any():
        print('Input arg "data_scrambled" should be a 2D numpy array of reals, with the same shape as "data_matrix"!')
        return 0
    if fft_axis not in [0, 1]:
        print('Input arg "fft_axis" should be 0 or 1!')
        return 0

    # if fft_axis != 0, transpose the data
    if fft_axis == 1:
        data_matrix = np.transpose(data_matrix)
        data_scrambled = np.transpose(data_scrambled)

    # init output list
    test_results = [False, False, False]

    # Check FFT component amplitudes / magnitudes. They should be the same, with differences only due to numeric inaccuracies
    # do forward FFT, use version for reals, treat data as if vars were in columns
    data_fft = np.fft.rfft(data_matrix, axis=0)
    data_scrambled_fft = np.fft.rfft(data_scrambled, axis=0)
    # compare magnitudes
    amp_diffs = np.abs(data_fft)-np.abs(data_scrambled_fft)
    print('Maximum difference between FFT component magnitudes: {:.3e}'.format(amp_diffs.max()))
    # set relevant output to True if passed the test
    if not (amp_diffs>epsilon).any():
        test_results[0] = True
        print('First test passed, original and scrambled data have matching FFT component amplitudes.')
    else:
        print('First test failed, found substantial difference\n' +
              'between original and scrambled data FFT component magnitudes.')

    # Check the covariance matrices
    data_cov = np.cov(data_matrix, rowvar=False)
    data_scrambled_cov = np.cov(data_scrambled, rowvar=False)
    maxDiff = (data_cov-data_scrambled_cov).flatten().max()  # maximum difference
    print('Maximum difference between covariance matrices: {:.3e}'.format(maxDiff))
    # set relevant output to True if passed the test
    if maxDiff <= epsilon:
        test_results[1] = True
        print('Second test passed, original and scrambled data have matching covariance structures.')
    else:
        print('Second test failed, found substantial difference\n' +
              'between original and scrambled data covariance structures.')

    # Check the correlations across original and scrambled vars
    ccoeffs = fastColumnCorr(data_matrix, data_scrambled)
    # if the mean and median are close to zero (<0.05) we consider that a success
    tmp = mplot.hist(ccoeffs, bins=data_matrix.shape[1]//40)
    if np.mean(ccoeffs) < 0.05 and np.median(ccoeffs) < 0.05:
        test_results[2] = True
        print('Third test passed, correlation coefficients group around 0.\n' +
              ' Look at the histogram for further details.')
    else:
        print('Third test failed, correlation coefficients seem to be biased.\n' +
              'Look at the histogram for further details.')

    return test_results