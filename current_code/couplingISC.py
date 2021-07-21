
import numpy as np
from numpy.random import default_rng


def random_normal_data(tr_no, voxel_no):
    """
    Helper function generating random standard normal data with given TR and voxel numbers.
    Output is a 2D numpy ndarray with shape (tr_no, voxel_no)
    """
    mu, sigma = 0, 1
    data = default_rng().normal(mu, sigma, (tr_no, voxel_no))
    return data


def timeshifted_data_1d(data, shift, padding_value='zero'):
    """
    Generates a matrix ("data_shifted") from an input vector ("data"),
    where the columns of the matrix are "shifted" versions of the input.
    The amount of shifting is from +shift" to -"shift" (that is, range(-shift, shift+1, 1)).

    Specifically, column 0 will contain the input vector data[shift:],
    padded with "padding_value" at the end, so that its length is
    the same as the length of "data".
    Column 1 then is data[shift-1:] plus padding, and so on,
    until the last column contains first padding, and data[0:-shift].


    Inputs
    data:           2D numpy array, with shape (n, 1), where n is the length of the input data.
    shift:          Positive integer. Maximum number of data points to shift.
                    E.g., if shift = 2, the columns of the output matrix will be
                    created from the input vector, shifted with range(-shift, shift+1, 1).
    padding_value:  String, either "zero" or "mean". Value for padding the input vector when it is shifted.
                    "mean" corresponds to the mean of the input vector.

    Output
    data_shifted:   2D numpy array with shape (n, shift*2+1), where n is the length of the input data.
                    Each column of "data_shifted" is generated from "data" by shifting it with a value from
                    range(-shift, shift+1, 1).
                    For example, with data = np.arange(10), shift = 2 and padding_value = 'zero',
                    data_shifted is:
                    array([[2., 1., 0., 0., 0.],
                           [3., 2., 1., 0., 0.],
                           [4., 3., 2., 1., 0.],
                           [5., 4., 3., 2., 1.],
                           [6., 5., 4., 3., 2.],
                           [7., 6., 5., 4., 3.],
                           [8., 7., 6., 5., 4.],
                           [9., 8., 7., 6., 5.],
                           [0., 9., 8., 7., 6.],
                           [0., 0., 9., 8., 7.]])

    """

    # check input data format
    if data.shape[1] != 1:
        raise ValueError('Input arg ''data'' should have shape (n, 1)!')
    # get padding value
    if padding_value == 'zero':
        pad = 0
    elif padding_value == 'mean':
        pad = np.mean(data)
    else:
        raise ValueError('Input arg ''padding_value'' should be either ''zero'' or ''mean''!')

    # preallocate output matrix
    data_shifted = np.zeros((data.shape[0], shift*2+1))

    # loop through data shifts
    for i in range(-shift, shift+1, 1):
        if i <= 0:
            data_shifted[:, i+shift] = np.concatenate((data[-i:, 0], np.tile(pad, -i)))
        else:
            data_shifted[:, i+shift] = np.concatenate((np.tile(pad, i), data[0:-i, 0]))

    return data_shifted


def timeshifted_data_2d(data, shift, padding_value='zero'):
    """
    Same as timeshifted_data_1d but with a 2D array as input. In other words, performs shifting of multiple data
    vectors (e.g. multiple voxel timeseries).

    Inputs
    data:           2D numpy array, with shape (n, v), where n is the length of a timeseries (e.g. TRs)
                    and v is the number of variables (e.g. voxels).
    shift:          Positive integer. Maximum number of data points to shift.
                    E.g., if shift = 2, the columns of the output matrix will be
                    created from the input vector, shifted with range(-shift, shift+1, 1).
    padding_value:  String, either "zero" or "mean". Value for padding the input vector when it is shifted.
                    "mean" corresponds to the mean of the input vector.

    Output
    data_shifted:   3D numpy array with shape (v, n, shift*2+1). See input arg "data" for dimensions.
                    Each column of "data_shifted" is generated from "data" by shifting it with a value from
                    range(-shift, shift+1, 1).
    """

    # check input data format
    if np.ndim(data) != 2:
        raise ValueError('Input arg ''data'' should be 2D!')
    # get padding value
    if padding_value == 'zero':
        pad = np.asarray([0])
    elif padding_value == 'mean':
        pad = np.mean(data, axis=0)  # returns vector
    else:
        raise ValueError('Input arg ''padding_value'' should be either ''zero'' or ''mean''!')

    # dimensions of data
    tr_no, vox_no = data.shape

    # preallocate output matrix
    data_shifted = np.zeros((vox_no, tr_no, shift*2+1))

    # loop through data shifts - different versions depending on the padding value

    # If padding is with zero (single value):
    if pad.shape[0] == 1:
        for i in range(-shift, shift+1, 1):
            if i <= 0:
                tmp = np.concatenate((data[-i:, :], np.tile(pad, (-i, vox_no))))
            else:
                tmp = np.concatenate((np.tile(pad, (i, vox_no)), data[0:-i, :]))
            data_shifted[:, :, i+shift] = tmp.T

    # If padding is with mean (vector):
    else:
        for i in range(-shift, shift+1, 1):
            if i <= 0:
                tmp = np.concatenate((data[-i:, :], np.tile(pad, (-i, 1))))
            else:
                tmp = np.concatenate((np.tile(pad, (i, 1)), data[0:-i, :]))
            data_shifted[:, :, i+shift] = tmp.T

    return data_shifted


def get_coupled_set_1d(tr_no, beta=None, shift=2, noise_level=0):
    """
    Helper function to generate an independent - dependent variable pair for testing coupling (OLS solution) methods.
    For given parameters, it generates a random independent data vector ("X", "speaker" data for our fMRI
    coupling use case) and a corresponding dependent data vector ("Y", "listener" data for our fMRI use case).
    "Y" is calculated by applying the "beta" linear coefficients to the time-shifted versions of "X".
    Additionally, (random standard normal) noise is added ("noise_level").

    Inputs
    tr_no:          Integer, length of data vectors (number of TRs in case of fMRI data).
    beta:           1D numpy array or list of values, with length "shift"*2+1.
                    Coefficients used for deriving the dependent variable from the independent.
                    Defaults to np.asarray([0.1, 0, 0.5, 2, 0]).
    shift:          Positive integer. Maximum number of data points to shift for the independent variable before
                    calculating the dependent variable.
                    E.g., if shift = 2, the independent data is shifted with range(-shift, shift+1, 1), then the
                    independent variable will be calculated as ("shifted_independent_var" @ "beta") (+ normalization).
    noise_level:    Positive number or zero. If not zero, a random standard normal vector ("noise") scaled by
                    "noise_level" is added to "Y" after it is calculated from "X" and "beta".

    Outputs:
    X:          1D numpy array. Independent variable data vector with shape (tr_no, 1)
    Y:          1D numpy array. Dependent variable data vector with shape (tr_no, 1)
    """

    # input checks
    if beta is None:
        beta = np.asarray([0.1, 0, 0.5, 2, 0])
    else:
        beta = np.asarray(beta)
        if np.ndim(beta) != 1:
            raise ValueError('Input arg ''beta'' should be a 1D array!')
    if shift % 1 != 0 or shift <= 0:
        raise ValueError('Input arg ''shift'' should be a positive integer!')
    if beta.shape[0] != shift*2+1:
        raise ValueError('Input arg ''beta'' should have length ''shift''*2+1!')
    if noise_level < 0:
        raise ValueError('Input arg ''noise_level'' should be a positive number or zero!')

    # generate data
    X = random_normal_data(tr_no, 1)
    X_shifted = timeshifted_data_1d(X, shift, padding_value='zero')
    Y = (X_shifted @ beta) / np.sum(beta)

    # add noise if requested
    if noise_level != 0:
        n = random_normal_data(tr_no, 1) * noise_level
        Y = (Y + n) / (1 + noise_level)

    return X, Y


def get_coupled_set_2d(data_dims, beta=None, shift=2, noise_level=0):
    """
    Same as get_random_set1d but generating multivariate sets (e.g. many voxels' worth of data
    at once for our fMRI use case). Outputs "X" and "Y" ("speaker" and "listener" data in our use case) are
    2D arrays with shape (tr_no, vox_no), otherwise it works the same way as the 1D version of the function.

    Inputs
    data_dims:       Tuple of integers, dimensions of data (number of TRs and voxels in case of fMRI data).
    beta:            1D numpy array or list of values, with length "shift"*2+1.
                     Coefficients used for deriving the dependent variable from the independent.
                     Defaults to np.asarray([0.1, 0, 0.5, 2, 0]).
    shift:           Positive integer. Maximum number of data points to shift for the independent variable before
                     calculating the dependent variable.
                     E.g., if shift = 2, the independent data is shifted with range(-shift, shift+1, 1), then the
                     independent variable will be calculated as ("shifted_independent_var" @ "beta") (+ normalization).
    noise_level:     Positive number or zero. If not zero, a random standard normal vector ("noise") scaled by
                     "noise_level" is added to "Y" after it is calculated from "X" and "beta".

    Outputs:
    X:              2D numpy array. Independent variable data vector with shape (tr_no, vox_no)
    Y:              2D numpy array. Dependent variable data vector with shape (tr_no, vox_no)
    """

    # input checks
    if beta is None:
        beta = np.asarray([0.1, 0, 0.5, 2, 0])
    else:
        beta = np.asarray(beta)
        if np.ndim(beta) != 1:
            raise ValueError('Input arg ''beta'' should be a 1D array!')
    if shift % 1 != 0 or shift <= 0:
        raise ValueError('Input arg ''shift'' should be a positive integer!')
    if beta.shape[0] != shift*2+1:
        raise ValueError('Input arg ''beta'' should have length ''shift''*2+1!')
    if noise_level < 0:
        raise ValueError('Input arg ''noise_level'' should be a positive number or zero!')
    try:
        tr_no, vox_no = data_dims
    except Exception:
        raise ValueError('Input arg ''data_dims'' should be a tuple of integers!')

    # generate data
    X = random_normal_data(tr_no, vox_no)
    X_shifted = timeshifted_data_2d(X, shift, padding_value='zero')
    Y = (X_shifted @ beta) / np.sum(beta)
    Y = Y.T

    # add noise if requested
    if noise_level != 0:
        n = random_normal_data(tr_no, 1) * noise_level
        Y = (Y + n) / (1 + noise_level)

    return X, Y


def coupling_loop(speaker, listener, shift, padding_value='zero'):
    """
    Coupling ISC estimation. The idea is to model the listener's BOLD time series as a
    linear function of the speaker's time shifted time series.
    See Stephens et al., 2010 and Silbert et al., 2014 for details.

    Written with fMRI data in mind.

    Inputs
    speaker:        2D numpy array, TRs X voxels (each column is a separate variable)
    listener:       2D numpy array, TRs X voxels (each column is a separate variable)
    shift:          Positive integer. Maximum number of data points to shift.
                    E.g., if shift = 2, the columns of the all speaker timeseries
                    will be shifted with the values in range(-shift, shift+1, 1).
    padding_value:  String, either "zero" or "mean". Value for padding speaker's data when it is shifted.
                    "mean" corresponds to the mean of each voxel timeseries.

    Outputs
    beta:           2D numpy array, voxels X shifts, the linear coefficients of shifted timeseries.
    residuals:      1D numpy array (column vector), voxels X 1. Sum of squared residuals
                    after solving the linear system in the least-squares sense.
    rsquared:       1D numpy array (column vector), voxels X 1. R squared values (goodness-of-fit measure)
                    for the coupling solution. Describes how well the speaker time series predict
                    the listener time series.
    """

    # input checks
    if speaker.shape != listener.shape or np.ndim(speaker) != 2 or np.ndim(listener) != 2:
        raise ValueError('Input args ''speaker'' and ''listener'' should have the same shape (as 2D arrays)!')
    if shift % 1 != 0 or shift <= 0:
        raise ValueError('Input arg ''shift'' should be a positive integer!')

    # get TR and voxel numbers
    tr_no, vox_no = speaker.shape

    # preallocate for coefficients (beta), sum of squared residuals,
    beta = np.zeros((vox_no, shift*2+1))
    residuals = np.zeros((vox_no, 1))
    rsquared = np.zeros((vox_no, 1))

    for i in range(vox_no):
        # for readability, define the data for given voxel
        speaker_vox = speaker[:, i]
        listener_vox = listener[:, i]
        # get time-shifted model (speaker) data for given voxel
        speaker_vox_shifted = timeshifted_data_1d(speaker_vox[:, np.newaxis], shift, padding_value)

        # Solve for coefficients, use standard least squares method.
        # An alternative way is to perform the calculation step-by-step ourselves:
        # X = speaker_vox_shifted
        # Y = listener_vox
        # b = (np.linalg.inv(X.T @ X) @ X.T) @ Y
        # We might be able to use this latter method on a stack of matrices,
        # that is, avoiding the for loop.
        ls_results = np.linalg.lstsq(speaker_vox_shifted, listener_vox, rcond=None)
        beta[i, :] = ls_results[0]
        residuals[i] = ls_results[1]
        # get R squared
        tmp_sstotal = np.sum(np.square((listener_vox - np.mean(listener_vox))))
        rsquared[i] = 1 - (ls_results[1]/tmp_sstotal)

    return beta, residuals, rsquared


def coupling_onestep(speaker, listener, shift, padding_value='zero'):
    """
    Same as the function 'coupling_loop' but with a hopefully faster implementation
    avoiding the for loop across voxels.
    The idea is to calculate the matrix multiplication and matrix inversion steps explicitly,
    as lower-level methods can be invoked on stacks of matrices directly.

    Inputs
    speaker:        2D numpy array, TRs X voxels (each column is a separate variable)
    listener:       2D numpy array, TRs X voxels (each column is a separate variable)
    shift:          Positive integer. Maximum number of data points to shift.
                    E.g., if shift = 2, the columns of the all speaker timeseries
                    will be shifted with the values in range(-shift, shift+1, 1).
    padding_value:  String, either "zero" or "mean". Value for padding speaker's data when it is shifted.
                    "mean" corresponds to the mean of each voxel timeseries.

    Outputs
    beta:           2D numpy array, voxels X shifts, the linear coefficients of shifted timeseries.
    residuals:      1D numpy array (column vector), voxels X 1. Sum of squared residuals
                    after solving the linear system in the least-squares sense.
    """

    # input checks
    if speaker.shape != listener.shape or np.ndim(speaker) != 2 or np.ndim(listener) != 2:
        raise ValueError('Input args ''speaker'' and ''listener'' should have the same shape (as 2D arrays)!')
    if shift % 1 != 0 or shift <= 0:
        raise ValueError('Input arg ''shift'' should be a positive integer!')

    # get shifted voxel-level time series for speaker data
    speaker_shifted = timeshifted_data_2d(speaker, shift, padding_value)  # returns 3D array, shape (v, tr, shifts)

    # get transposed matrices for each voxel
    speaker_shifted_T = np.transpose(speaker_shifted, axes=[0, 2, 1])  # 3D array, shape (v, shifts, tr)
    listener_T = np.transpose(listener, axes=[1, 0])[:, :, np.newaxis]  # 3D array, shape (v, tr, 1)

    # the least-squares solution is given by: b = (X.T @ X)^(-1) @ X.T @ Y
    beta = (np.linalg.inv(speaker_shifted_T @ speaker_shifted) @ speaker_shifted_T) @ listener_T

    # get residuals


    return np.squeeze(beta)


def coupling_test(data_dims, shift=2, coupling_function='loop'):
    """
    Test the coupling method. Main steps:
    (1) Generate random coupled data set with known coefficients, using get_coupled_set_2d
    (2) Estimate coupling with 'coupling_loop' or 'coupling_onestep'
    (3) Compare the results to the known coefficients

    Inputs
    data_dims:          Tuple of integers, dimensions of both speaker and listener data
                        (number of TRs and voxels in case of fMRI data).
    shift:              Positive integer. Maximum number of data points to shift.
                        E.g., if shift = 2, the columns of the all speaker timeseries
                        will be shifted with the values in range(-shift, shift+1, 1).
    coupling_function:  String, either 'loop' or 'onestep'. Defines the coupling estimation function to use.

    Output
    beta_known:         1D numpy array, the known coefficients
    beta_est:           1D numpy array, estimated coefficients
    """

    # input checks
    if shift % 1 != 0 or shift <= 0:
        raise ValueError('Input arg ''shift'' should be a positive integer!')
    try:
        tr_no, vox_no = data_dims
    except Exception:
        raise ValueError('Input arg ''data_dims'' should be a tuple of integers!')

    # generate beta vector
    beta_known = default_rng().random((shift*2+1, 1))*2-1  # values in range -1 : +1
    beta_known = beta_known/np.sum(beta_known)
    beta_known = beta_known[:, 0]

    # generate random coupled set
    speaker, listener = get_coupled_set_2d(data_dims, beta=beta_known, shift=shift, noise_level=0)

    # estimate coupling
    if coupling_function == 'loop':
        beta_est, residuals, rsquared = coupling_loop(speaker, listener, shift)
    elif coupling_function == 'onestep':
        beta_est = coupling_onestep(speaker, listener, shift)

    # compare coupling
    beta_equal = np.zeros((vox_no, 1))
    for i in range(vox_no):
        beta_equal[i, 0] = np.allclose(beta_known, beta_est[i, :])

    # report results
    if np.sum(beta_equal) == vox_no:
        print('COEFFS ARE EQUAL, TEST PASSED!')
    else:
        print('COEFFS ARE NOT EQUAL, TEST FAILED!')

    return beta_known, beta_est
