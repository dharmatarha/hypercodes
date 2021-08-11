import numpy as np
import statsmodels.stats.multitest as multi

def coupling_lm_permTest_v2(realData, permData, alpha):

    """
    This function is a supplement to the coupling lag model scripts written for Thalia's hyperscanning dataset.
    It takes a vector of real model parameters and uses a 2d array of parameters estimated after circle shifting
    speaker time series to generate a null distribution and get p-values. It also extracts some summary measures
    and mapping variables to facilitate downstream stat map generation (see output description below).

    :param realData:           2D numpy array of real model parameters [voxels x 1]
    :param permData:           2D numpy array [permutations x voxels]
    :param alpha:              alpha level at which to evaluate hypothesis tests
    :return:                   structures containing several FDR-related outputs:
             pVal:             p-values
                   right_tail: proportion of null distribution to the RIGHT of the real parameter estimate [voxels x 1]
                   left_tail:  proportion of null distribution to the LEFT of the real parameter estimate [voxels x 1]
                   two_tail:   minimum value between pVal_right_tail and pVal_left_tail [voxels x 1]
                   map:        boolean vector to indicate the sidedness of pVal_two_tail values -- True = left-tailed
                               p-value, False = right-tailed p-value
                   fdr_right:  FDR corrected version of pVal_right_tail -- two vectors with length = # voxels, [0]:
                               boolean where True means corrected p-value less than alpha, [1]: corrected p-values
                   fdr_left:   FDR corrected version of pVal_left_tail
                   fdr_right:  FDR corrected version of pVal_two_tail (evaluated at alpha / 2)
             mask:             realData filtered such that voxels that don't survive FDR correction are set to zero.
                   right:      right-tailed test
                   left:       left-tailed test
                   two:        two-tailed test
                   two_above:  two-tailed test filtering out any significant hits on the left side of the null
                               (equivalent to a right-tailed test at alpha / 2)
                   two_below:  two-tailed test filtering out any significant hits on the left side of the null
                               (equivalent to a left-tailed test at alpha / 2)
             pSig:             proportions of voxels with FDR-corrected significant p-values
                   right:      right-tailed test
                   left:       left-tailed test
                   two:        two-tailed test
                   two_above:  proportion of two-tailed test significant voxels on the right side of null distribution
                   two_below:  proportion of two-tailed test significant voxels on the left side of null distribution


    """

    # check inputs
    if not isinstance(realData, np.ndarray) or realData.shape[1] != 1:
        print('Input arg "realData" should be a 2D numpy array of reals with dimensions of [voxels x 1]!')
        if len(realData.shape) == 1:
            print('We have detected that "realData" was entered as a 1D vector -- a totally reasonable mistake.')
            print('Adding a second dimension...')
            realData = realData.reshape([len(realData),1])
        return 0

    if not isinstance(permData, np.ndarray) or len(permData.shape) != 2:
        print('Input arg "permData" should be a 2D numpy array of reals!')
        return 0

    if not len(realData) == permData.shape[0]:
        print('The length of "realData" should be equal to the number of rows in "permData"!')
        if len(realData) == permData.shape[1]:
            print('We have detected that the number of columns in "permData" is equal to the length of "realData".')
            print('Assuming that "permData" was inputted with permutations as columns... transposing the "permData" matrix...')
            permData = np.transpose(permData)
        else:
            return 0

    # get number of voxels and permutations from input arrays
    nVox = len(realData)
    nPermut = permData.shape[1]

    # initialize output dictionaries
    pVals = {'right_tail':np.empty(nVox),
             'left_tail':np.empty(nVox),
             'two_tail':np.empty(nVox),
             'map':np.empty(nVox),
             'fdr_right':np.empty(nVox),
             'fdr_left':np.empty(nVox),
             'fdr_two':np.empty(nVox)}
    mask = {'right':np.empty(nVox),
            'left':np.empty(nVox),
            'two':np.empty(nVox),
            'two_right':np.empty(nVox),
            'two_left':np.empty(nVox),}
    pSig = {'right': np.empty(nVox),
            'left': np.empty(nVox),
            'twp': np.empty(nVox),
            'two_above': np.empty(nVox),
            'two_below': np.empty(nVox)}

    # get permutation test p-values
    pVals['right_tail'] = np.sum((permData > realData), axis=1) / nPermut
    pVals['left_tail'] = np.sum((permData < realData), axis=1) / nPermut
    pVals['two_tail'] = np.minimum(pVals['right_tail'], pVals['left_tail'])
    pVals['map'] = pVals['left_tail'] > pVals['right_tail']

    # FDR correction - outputs two vectors of length = # voxels, [0]: boolean where True means corrected p-value
    # less than alpha, [1]: corrected p-values
    pVals['fdr_right'] = multi.fdrcorrection(pVals['right_tail'], alpha=alpha)
    pVals['fdr_left'] = multi.fdrcorrection(pVals['left_tail'], alpha=alpha)
    pVals['fdr_two'] = multi.fdrcorrection(pVals['two_tail'], alpha=alpha / 2)

    # get R^2 masks where all voxels that did not survive FDR correction are set to zero
    mask['right'] = np.copy(realData)
    mask['right'][np.invert(pVals['fdr_right'][0])] = 0
    mask['left'] = np.copy(realData)
    mask['left'][np.invert(pVals['fdr_left'][0])] = 0
    mask['two'] = np.copy(realData)
    mask['two'][np.invert(pVals['fdr_two'][0])] = 0
    mask['two_above'] = np.copy(realData)
    mask['two_above'][np.minimum(pVals['fdr_two'][0], np.invert(pVals['map']))] = 0
    mask['two_below'] = np.copy(realData)
    mask['two_below'][np.minimum(pVals['fdr_two'][0], pVals['map'])] = 0

    # proportion of voxels that show significant FDR CORRECTED parameter estimates
    pSig['right'] = len(np.where(pVals['fdr_right'][0])[0]) / nVox # right-tailed
    pSig['left'] = len(np.where(pVals['fdr_left'][0])[0]) / nVox # left-tailed
    pSig['two'] = len(np.where(pVals['fdr_two'][0])[0]) / nVox # any two-tailed
    pSig['two_above'] = np.count_nonzero(np.minimum(pVals['fdr_two'][0], np.invert(pVals['map']))) / nVox  # two-tailed on right side
    pSig['two_below'] = np.count_nonzero(np.minimum(pVals['fdr_two'][0], pVals['map'])) / nVox # two-tailed on left side

    return pVals, mask, pSig