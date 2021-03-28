
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


def get_random_set(tr_no):
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
            data_shifted[:, i+shift] = np.concatenate((data[-i:], np.tile(pad, -i)))
        else:
            data_shifted[:, i+shift] = np.concatenate((np.tile(pad, i), data[0:-i]))

    return data_shifted


def timeshifted_data_2d(data, shift, padding_value ='zero'):
    """
    Same as timeshifted_data_1d but with a 2D array as input. In other words, performs shifting of multiple data
    vectors (e.g. multiple voxel timeseries).
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

