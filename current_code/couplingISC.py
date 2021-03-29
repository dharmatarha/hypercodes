
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


def get_random_set_1d(tr_no):
    """
    Helper function to generate an independent - dependent variable pair for testing OLS solution methods
    """
    b = np.asarray([0.1, 0, 0.5, 2, 0])
    tmp = random_normal_data(tr_no, 1)
    X = timeshifted_data_1d(tmp, 2, padding_value='zero')
    Y = np.matmul(X, b)
    return X, b, Y



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
    data:           1D numpy array, with shape (n, 1), where n is the length of the input data.
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
    if data.shape[1] != 1 :
        raise ValueError('Input arg ''data'' should have shape (n, 1)!')
    # get padding value
    if padding_value is 'zero':
        pad = 0
    elif padding_value is 'mean':
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


def timeshifted_data_2d(data, shift, padding_value ='zero'):
    """
    Same as timeshifted_data_1d but with a 2D array as input. In other words, performs shifting of multiple data
    vectors (e.g. multiple voxel timeseries).

    Inputs
    data:           2D numpy array, with shape (n, v), where n is the length of a timeseries
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

    return data_shifted


def coupling(speaker, listener, shift, padding_value='zero'):
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
    """

    # simple version, looping through voxels:
    tr_no, vox_no = speaker.shape

    # preallocate coeffs, residuals,
    betas = np.zeros((shift*2+1, vox_no))
    residuals = np.zeros((vox_no, 1))

    for i in range(vox_no):
        # for readability, define the data for given voxel
        speaker_vox = speaker[:, i]
        listener_vox = listener[:, i]
        # get time-shifted model (speaker) data for given voxel
        speaker_vox_shifted = timeshifted_data_1d(speaker_vox, shift, padding_value)

        # Solve for coefficients, use standard least squares method.
        # An alternative way is to perform the calculation step-by-step ourselves:
        # X = speaker_vox_shifted
        # Y = listener_vox
        # b = (np.linalg.inv(X.T @ X) @ X.T) @ Y
        # We might be able to use this latter method on a stack of matrices,
        # that is, avoiding the for loop.

        ls_results = np.linalg.lstsq(speaker_vox_shifted, listener_vox, rcond=None)
        betas[:, i] = ls_results[0]
        residuals[i] = ls_results[1]

    return betas, residuals