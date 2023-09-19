import numpy as np
from scipy.spatial.distance import cdist

def DIANA(data, k):
    clustAssing = []
    n_objects = data.shape[0]
    clusters = [list(range(n_objects))]
    # dists = distance_matrix(objects, objects, p=2)
    objects = data
    # transpose
    at = objects.T
    # covariance matrix
    D = np.cov(at)
    invD = np.linalg.inv(D)
    dists = cdist(objects, objects, metric='mahalanobis',VI=invD)
    # print(type(dists))
    while len(clusters) < k:
        print("Length of clusters = {}".format(len(clusters)))
        diameters = [np.max(dists[ctr][:, ctr]) for ctr in clusters]
        max_dtr = np.argmax(diameters)

        avg_dists_in_ctr = np.mean(dists[clusters[max_dtr]][:, clusters[max_dtr]], axis=1)
        seed_idx = np.argmax(avg_dists_in_ctr)

        splinter_group = [clusters[max_dtr][seed_idx]]
        non_splinter = clusters[max_dtr].copy()
        non_splinter.pop(seed_idx)

        while True:
            in_dist = np.mean(dists[non_splinter][:, splinter_group], axis=1)
            out_dist = dists[non_splinter][:, non_splinter]
            # to exclude i to i distance(=1)
            out_dist = np.mean(out_dist, axis=1) - 1 / len(non_splinter)
            d_h = out_dist - in_dist
            idx = np.argmax(d_h)
            if d_h[idx] > 0:
                splinter_group.append(non_splinter[idx])
                non_splinter.pop(idx)
            else:
                break
        del clusters[max_dtr]
        clusters.append(non_splinter)
        clusters.append(splinter_group)

    for i in range(k):
        clustAssing.append(clusters[i])
    return clustAssing