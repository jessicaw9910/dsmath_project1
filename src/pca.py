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

def plot_pca(eigenvalues, n_pc, title='pca', folder='plots/', save=True, img_size=(12,4)):
    fig = plt.figure(figsize=img_size);
    ## variance
    plt.subplot(1, 2, 1);
    plt.plot(range(len(eigenvalues)), eigenvalues);
    plt.vlines(n_pc, ymin=0, ymax=np.max(eigenvalues), colors='r', linewidths=2, linestyle='dashed')
    plt.xlabel('$\it{k}$th PC', fontsize=12);
    plt.ylabel('$\it{\sigma_k^2}$', fontsize=12);
    ## percent variance explained
    plt.subplot(1, 2, 2);
    plt.plot(range(len(eigenvalues)), np.cumsum(eigenvalues / np.sum(eigenvalues)));
    plt.vlines(n_pc, ymin=0, ymax=1, colors='r', linewidths=2, linestyle='dashed')
    plt.xlabel('$\it{k}$th PC', fontsize=12);
    plt.ylabel('% $\it{\sigma_k^2}$ explained by $\it{k}$th PC', fontsize=12);
    if save:
        plt.savefig(folder + title + '.svg',bbox_inches='tight');
    plt.show();