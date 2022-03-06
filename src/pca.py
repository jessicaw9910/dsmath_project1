#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def find_pc(mx_input):
    ''' Find the principal component values and directions
    
    Params:
    - mx_input: input matrix

    Output:
    - eigenvalues_sorted: 
    - eigenvectors_sorted:
    '''

    # calculate covariance of mean-centered data
    covar_mx = np.cov(mx_input, rowvar=False)

    # used eigh instead of eig since symmetric
    eigenvalues, eigenvectors = np.linalg.eigh(covar_mx)

    # sort decreasing order
    idx_sort = np.argsort(eigenvalues)[::-1]
    eigenvalues_sorted = eigenvalues[idx_sort]
    eigenvectors_sorted = eigenvectors[:, idx_sort]
    
    return eigenvalues_sorted, eigenvectors_sorted

def project_pca(mx_input, eigenvectors, n_pc=10):
    eigenvectors_subset = eigenvectors[:, 0:n_pc]
    mx_pca = (eigenvectors_subset.T @ mx_input.T).T
    return mx_pca

def plot_pca(eigenvalues, img_size=(12,4), title='pca', folder='plot/', save=True):
    fig = plt.figure(figsize=img_size);
    ## variance
    plt.subplot(1, 2, 1);
    plt.plot(range(len(eigenvalues)), eigenvalues);
    plt.xlabel('$\it{k}$th principal component');
    plt.ylabel('$\it{\sigma_k^2}$');
    ## percent variance explained
    plt.subplot(1, 2, 2);
    plt.plot(range(len(eigenvalues)), np.cumsum(eigenvalues / np.sum(eigenvalues)));
    plt.xlabel('$\it{k}$th principal component');
    plt.ylabel('% $\it{\sigma_k^2}$ explained by $\it{k}$th PC');
    if save:
        plt.savefig(folder + title + '.pdf',bbox_inches='tight');
    plt.show();