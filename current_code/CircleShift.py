#necessary imports
import numpy as np
import random


"""
Script documentation.

Set of functions that can be used to circle shift data and test its performance.

Functions:
-Circle_Shift: circlularly shift a 2D np array to generate permuted timeseries
-array2autocorrelation: find the lag-n autocorrelation in each variable in a 2d np array 
-autocorrelation_test: compare if autocorrelation present in two 2d arrays of variables 
 is similar up to a set thershold
-generate_ac_data: generate a 2d array of data with predefined amounts of autocorrelaiton 
 per variable (column)

Other good tests for permuted data can be found in Adam's phaseScramble script (phase_scrambling_tests)

"""


def Circle_Shift(data_matrix, shift_point = float("nan"), shift_axis = 1, rand_shift = True):
    """
    Circle shifts an input data matrix. Essentially takes the bottom n rows of a matrix and tosses them at the top.
    The variable n is determined by the shift_point variable (or, to randomly select it with each call,
    set rand_shift to true). 
    
    This approach is a) suggested for use in ISC analyses by Nastase et al., 2019, 
    b) used in an Intersubject Correlations tutorial by Juha Lahnakoski and Luke Chang, and
    c) described in section 7.4 of Lancaster et al., 2008.
    
    Input requirements:
    -data_matrix: 2D np array with rows = samples and columns = variables.
    -shift_point: int between 0 and the max number of samples (only matters if rand_shift = False)
    -shift_axis: either 0 or 1 (ints)
    -rand_shift: Boolean (True or False)
    -numpy loaded as np
    -random loaded as random
    
    Input details:
    -If rand_shift is set to true, shift points will be automatically picked at random with each function call.
     Points are picked within the range of (0,number of samples) with the endpoints excluded
    -If rand_shift is set to false, you must supply the shift point via the shift point variable. 
    -The matrix is automatically shifted along columns (i.e. the bottom n rows are moved to the top). 
     If your data has rows = variables and columns = samples, then set shift_axis = 0.
    
    Returns:
    -A shifted version of data_matrix with the bottom n samples put at the start
    """
    
    ### Preprocessing and input checks ###
    
    #check those inputs
    if not isinstance(data_matrix, np.ndarray) or len(data_matrix.shape) != 2:
        print('Input arg "data_matrix" should be a 2D numpy array of reals!')
        return 0
    if shift_axis not in [0, 1]:
        print('Input arg "shift_axis" should be 0 or 1!')
        return 0
    if rand_shift not in [True, False]:
        print('Input arg "rand_shift" should be True or False!')
        return 0
    
    #transpose data_matrix if shift_axis = 0 -> rows will now be samples, cols = vars
    if shift_axis == 0:
        data_matrix = np.transpose(data_matrix)
    
    #select a random shift value in the range (0,nsamples) if rand_shift = True. 
    if rand_shift == True:
        shift_point = random.randint(1,np.size(data_matrix,0)-1)
        
    #If rand_shift = False, make sure shift_point is an int in the specified range
    elif isinstance(shift_point, int) == False or shift_point <= 0 or shift_point >= np.size(data_matrix,0):
        print('Input arg "shift_point" must be an int in the range of (0,number of samples) - endpoints excluded')
        return 0
    
    ###The meat of the function###
    
    #take the bottom n rows and move them to the top
    shifted_data_matrix = np.vstack((data_matrix[shift_point:],data_matrix[0:shift_point]))
    
    #return data in desired format if axis switched
    if shift_axis == 0:
        shifted_data_matrix = np.transpose(shifted_data_matrix)
    
    return shifted_data_matrix


def array2autocorrelation(data_matrix, shift_axis=1, ac_lag=1):
    """
    Computes autocorrelaiton for each variable in data_matrix (2d numpy array assumed to be vars X samples).
    The lag at which autocorrelation is calculated is determined by ac_lag (default = 1).

    Inputs:
    -data_matrix: 2D np array with rows = samples and columns = variables. 
    -shift_axis: either 0 or 1 (ints)
    -ac_lag: must be an int greater than 0. Specifies the lag at which autocorrelation is calculated

    Input details:
    -data is assumed to be samples X variables. If your data has the shape variables X samples, set 
    shift_axis to 0.

    Returns:
    -one dimensional numpy array containing autocorrelation values across variables
    """

    #check those inputs
    if not isinstance(data_matrix, np.ndarray) or len(data_matrix.shape) != 2:
        print('Input arg "data_matrix" should be a 2D numpy array of reals!')
        return 0
    if shift_axis not in [0, 1]:
        print('Input arg "shift_axis" should be 0 or 1!')
        return 0
    if not isinstance(ac_lag, int) or ac_lag <= 0:
        print('Input arg "ac_lag" must be an int greater than 0!')
        return 0

    #transpose data_matrix if shift_axis = 0 -> rows will now be samples, cols = vars
    if shift_axis == 0:
        data_matrix = np.transpose(data_matrix)


    #get lagged autocorrelation for each variable in the real data (lag specified by ac_lag)

    #### USING LIST COMP FOR THIS - COULD BE MADE FASTER ####

    #steps through columns (vars) and finds correlation between start:end-lag and start+lag:end vectors
    #returns as a np.array
    return np.array([np.corrcoef(data_matrix[:-ac_lag,col], data_matrix[ac_lag:,col])[0][1] for col in range(data_matrix.shape[1])])

def autocorrelation_test(data_matrix, data_scrambled, shift_axis=1, max_diff= 0.01):
    """
    Compares observed lag-1 autocorrelaiton in natural data to lag-1 autocorrelation in scrambled data. 
    If the maximum difference in autocorrelation values across variables is less than a user-specified 
    criterion (max_diff = 0.01 as default), then the test passes.

    Inputs:
    -data_matrix: 2D np array with rows = samples and columns = variables.
    -data_scrambled: comparison data of the same shape as data_matrix. Ideally should be a scrambled or
     permuted version of data_matrix. 
    -shift_axis: either 0 or 1 (ints)
    -max_diff: criterion to determine if test passes (minimum acceptable difference in autocorrelation between 
    real and shuffled data). Should be a float in the range of 0 < max_diff < 1

    Input details:
    -data is assumed to be samples X variables. If your data has the shape variables X samples, set 
    shift_axis to 0.

    Returns:
    -one dimensional numpy array containing the difference in autocorrelation values across variables
    """

    ### Preprocessing and input checks ###
    
    #check those inputs
    if not isinstance(data_matrix, np.ndarray) or len(data_matrix.shape) != 2:
        print('Input arg "data_matrix" should be a 2D numpy array of reals!')
        return 0
    if not isinstance(data_scrambled, np.ndarray) or len(data_scrambled.shape) != 2:
        print('Input arg "data_scrambled" should be a 2D numpy array of reals!')
        return 0
    if data_matrix.shape != data_scrambled.shape:
        print('Real and scrambled data array sizes must match!')
        return 0 
    if shift_axis not in [0, 1]:
        print('Input arg "shift_axis" should be 0 or 1!')
        return 0
    if not isinstance(max_diff, float):
        print('Input arg max_diff must be a float')
        return 0 
    elif max_diff <= 0 or max_diff >=1:
        print('Input arg max_diff must be in the range of 0 < max_diff < 1')
        return 0
    
    #transpose data_matrix if shift_axis = 0 -> rows will now be samples, cols = vars
    if shift_axis == 0:
        data_matrix = np.transpose(data_matrix)

    #use array2autocorrelation to get ac values for real and scrambled data and finds abs of diffs between them
    ac_diffs = np.absolute(array2autocorrelation(data_matrix) - array2autocorrelation(data_scrambled))

    #determines if max difference is less than max_diff
    if np.max(ac_diffs) >= max_diff:
        print(f"Test failed! Maximum observed difference of {np.max(ac_diffs)} >= criterion of {max_diff}!")
    else:
        print(f"Test passed! Maximum observed difference of {np.max(ac_diffs)} < criterion of {max_diff}!")

    #return the np array of differences 
    return ac_diffs

def generate_ac_data(nVars, nSamples, autoCorr, mu=0, sigma=1):
    """
    Generates a random dataset with prespecified amounts of lag-1 autocorrelation per variable.

    Borrows code from this stack overflow response:
    https://stackoverflow.com/questions/33898665/python-generate-array-of-specific-autocorrelation

    Inputs:
    -nVars: int. the number of random variables (columns) to be generated
    -nSamples: int. the number of random samples (rows) from each variable to be generated
    -autoCorr: float or 1d numpy array of floats in the same range. If entering an array, the array
     must have one entry for each random variable (i.e. len(autoCorr) = nVars). All floats must be 
     the range -1 < autoCorr < 1. If only a single float is used, it is assumed that each variable has
     the same lag-1 autocorrelation.
    -mu: float, int or 1d numpy array of floats. The mean of each variable. 
     If only a single float or int is used, it is assumed that each variable has the same mean. 
     If entering an array, the array must have one entry for each random variable (i.e. len(mu) = nVars).
    -sigma: float, int or 1d numpy array of floats. The variance of each variable. 
     If only a single float or int is used, it is assumed that each variable has the same variance. 
     If entering an array, the array must have one entry for each random variable (i.e. len(sigma) = nVars).

    Returns:
    -2d numpy array of the shape samples X variables 
    """
    #check inputs!
    if not isinstance(nVars, int):
        print("nVars must be an int!")

    if not isinstance(nSamples, int):
        print("nSamples must be an int!")

    if isinstance(mu, int) or isinstance(mu, float):
        mu_float = mu
        #create array with same mu value for each entry (one per variable)
        mu = np.full(shape=nVars, fill_value=mu_float, dtype=np.float)
    elif not isinstance(mu, np.ndarray) or len(mu.shape) != 1 or len(mu) != nVars:
        print("mu must be a float or a 1d numpy array with one entry for each variable")
        return(0)

    if isinstance(sigma, int) or isinstance(sigma, float):
        sigma_float = sigma
        #create array with same sigma value for each entry (one per variable)
        sigma = np.full(shape=nVars, fill_value=sigma_float, dtype=np.float)
    elif not isinstance(sigma, np.ndarray) or len(sigma.shape) != 1 or len(sigma) != nVars:
        print("sigma must be a float or a 1d numpy array with one entry for each variable")
        return(0)

    if isinstance(autoCorr, float):
        autoCorr_float = autoCorr
        #create array with same autoCorr value for each entry (one per variable)
        autoCorr = np.full(shape=nVars, fill_value=autoCorr_float, dtype=np.float)
    elif not isinstance(autoCorr, np.ndarray) or len(autoCorr.shape) != 1 or len(autoCorr) != nVars:
        print("autoCorr must be a float or a 1d numpy array with one entry for each variable")
        return(0)

    #Check for ac values outside of specified range
    for i in autoCorr:
        if i <= -1 or i >= 1:
            print(f"Found an autoCorr of value: {i}.\nAll values of autoCorr must be in the range -1 < autocorr < 1")
            return(0)
    
    #If code reaches here, checks have passed. Time to make data.

    #Overview: Create list of variables. Convert to array. Transpose and return.

    #Loop through nVars and create 1d np array for each

    #storage list
    data_list = []

    for idx in range(nVars):
        
        #START OF MOSTLY BORROWED CODE

        # Find out the offset `c` and the std of the white noise `sigma_e`
        # that produce a signal with the desired mean and variance.
        # See https://en.wikipedia.org/wiki/Autoregressive_model
        # under section "Example: An AR(1) process".
        c = mu[idx] * (1 - autoCorr[idx])
        sigma_e = np.sqrt((sigma[idx] ** 2) * (1 - autoCorr[idx] ** 2))

        #Sample the auto-regressive process.
        signal = [c + np.random.normal(0, sigma_e)]
        for _ in range(1, nSamples):
            signal.append(c + autoCorr[idx] * signal[-1] + np.random.normal(0, sigma_e))

        #END OF MOSTLY BORROWED CODE

        data_list.append(np.array(signal))
    
    #return data_list as a np array and transpose it so rows = samples and cols = vars
    return np.transpose(np.array(data_list))
