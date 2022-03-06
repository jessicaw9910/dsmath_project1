#!/usr/bin/env python3

import numpy as np

def fit_svd(matrix_init, index_train, rank=5, n_iter=100, thresh=1e-4, verbose=False):
    '''Fit SVD model iteratively

    Params:
    matrix_init: imputed, zero-mean initialization
    index_train: index of training data
    rank: the rank of low rank model
    n_iter: the max number of iterations
    thresh: the threshold of convergence to 0
    verbose:   
    '''
    ## initialize yhat 
    yhat = matrix_init
    
    ## iterate over n_iter
    for ind in range(n_iter):

        ## SVD and then calc low-rank estimate
        U, s, V = np.linalg.svd(yhat, full_matrices=False)
        yhat = U[:,:rank] @ np.diag(s)[:rank,:rank] @ V[:rank,:]
        
        ## calculate RMSE
        y_train = matrix_init.flat[index_train]
        yhat_train = yhat.flat[index_train]
        RMSE = np.sqrt(np.mean((y_train - yhat_train)**2))
        
        ## interim results if verbose
        if verbose:
            print("Iteration " + str(ind) + " Error: " + str(RMSE))
        
        ## stopping criteria
        if (RMSE < thresh):
            break
            
    return RMSE, yhat, U, s, V