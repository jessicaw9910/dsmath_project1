#!/usr/bin/env python3

from sklearn.cluster import KMeans
from kneed import KneeLocator
from sklearn.metrics import silhouette_score

def find_kmeans(mx_input, dict_kwarg, n_clust=30):
    sse = []
    silhouette = []

    for k in range(1, n_clust+1):
        kmeans = KMeans(n_clusters=k, **dict_kwarg)
        kmeans.fit(mx_input)
        sse.append(kmeans.inertia_)
        if k > 1:
            score = silhouette_score(mx_input, kmeans.labels_)
            silhouette.append(score)

    kl = KneeLocator(range(1, n_clust+1), sse, curve="convex", direction="decreasing")
    n_clust = kl.elbow

    kmeans = KMeans(n_clusters=n_clust, **dict_kwarg)
    kmeans.fit(mx_input)
    
    return kmeans, sse, silhouette