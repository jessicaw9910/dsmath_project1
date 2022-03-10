# Analysis of GDSC Data

A brief video summary of this project can be found [here](https://www.youtube.com/watch?v=2u5BWscHOeg) on YouTube. 

## /assets

Contains relevant csv files

### Raw data

+ **GDSC2_fitted_dose_response_25Feb20.csv** - contains the raw IC<sub>50</sub> data for 135,242 drug / cell line combinations; data can be downloaded [here](ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/GDSC2_fitted_dose_response_25Feb20.xlsx)
+ **cell_list.csv** - contains GDSC-generated information about the included cell lines, including tissue and TCGA classification; data can be downloaded [here](ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/Cell_Lines_Details.xlsx)
+ **drug_data.csv** - contains GDSC-generated information about the included drugs, including drug pathways and targets; data can be found [here](https://www.cancerrxgene.org/downloads/drug_data)

### Processed data

+ **cell.csv** - results of the analysis from `GDSC_Project.ipynb` for the cell data, including contains two dimensions from PCA and t-SNE for plotting, the cluster identities from *k*-nearest neighbors using full, PCA-transformed, and low rank approximations of the data, and the mean lnIC<sub>50</sub> per cell line
+ **cell_lrm.csv** - low-rank approximation of the cell matrix

![formula](https://render.githubusercontent.com/render/math?math=\color{white}\large\U \Sigma^{\frac{1}{2}})

+ **drug.csv** - results of the analysis from `GDSC_Project.ipynb` for the drug data, including contains two dimensions from PCA and t-SNE for plotting, the cluster identities from *k*-nearest neighbors using full, PCA-transformed, and low rank approximations of the data, and the mean lnIC<sub>50</sub> per compound
+ **drug_lrm.csv** - low-rank approximation of the drug matrix

![formula](https://render.githubusercontent.com/render/math?math=\color{white}\large\Sigma^{\frac{1}{2}} V^T)

## /src

Contains relevant scripts and notebooks

* **GDSC_Project.R** - contains
* **GDSC_Project.ipynb** - contains 
* **kmeans.py**
    * `find_kmeans`: find an optimal number clusters via the elbow method and fit *k*-means with this many clusters
    * `plot_kmeans`: plot the SSE vs. clusters and elbow point for cell and drug data
* **lowrank.py**
    * `fit_svd`: fit SVD model iteratively for a given rank *r*
* **pca.py**
    * `find_pc`: returns the eigenvalues and eigenvectors of the covariance matrix
    * `project_pca`: transforms the original matrix via projection using a specified number of principal components
    * `plot_pca`: plot the variance by principal component and cumulativev variance by principal component
* **utils.py**
    * `import_data`: loads GDSC data and pre-process into a wide matrix
    * `process_data`: produces mean-centered data and masks

## Summary

![formula](https://render.githubusercontent.com/render/math?math=\color{white}\large\D_{i, j} \approx \sum_{l=1}^ra_l\[i\]b_l\[j\] = U \Sigma V^T)
![formula](https://render.githubusercontent.com/render/math?math=\color{white}\large\U \Sigma^{\frac{1}{2}})
![formula](https://render.githubusercontent.com/render/math?math=\color{white}\large\Sigma^{\frac{1}{2}} V^T)