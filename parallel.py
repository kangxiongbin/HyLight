import sys
import os
from multiprocessing import Pool

from assembly import get_superead4


def run_on_local2(num_clusters, cluster_dir, max_tip_len, rm_trans, threads):
    params = []
    for j in range(num_clusters):
        i = str(j + 1)
        param = [i, cluster_dir, max_tip_len, rm_trans]
        params.append(param)

    pool = Pool(threads)
    pool.map(get_superead4, params, chunksize=1)  # ordered
    pool.close()
    pool.join()
    return

def run_on_local3(num_clusters, cluster_dir, max_tip_len, rm_trans, threads):
    params = []
    for j in range(num_clusters):
        i = str(j + 1)
        param = [i, cluster_dir, max_tip_len, rm_trans]
        params.append(param)

    with Pool(threads) as pool:
        pool.map(get_superead4, params, chunksize=1)
    return