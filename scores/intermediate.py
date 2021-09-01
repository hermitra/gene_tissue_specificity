import numpy as np


# intermediate functions

def dpm(x):
    # dispersion measure of a vector
    n = len(x)
    if n == 1:
        return 0
    else:
        disp_msr = np.std(x, ddof=1) * np.sqrt(n)
        return disp_msr


def js_distance(p, q):
    # jensen-shannon distance (metric)
    p = p / np.sum(p)
    q = q / np.sum(q)

    left = entropy((p + q) / 2)
    right = (entropy(p) + entropy(q)) / 2

    js = left - right
    jsd = np.sqrt(js)
    return jsd


def entropy(x):
    # entropy
    n = len(x)
    if not np.any(x):
        return np.log2(n)
    else:
        p = x / np.sum(x)
        p = p[p != 0]
        h = -1 * np.dot(p, np.log2(p))
        return h


def tukey_biweight(x, c=5, epsilon=1e-4):
    # measure of central tendency: tukey biweight
    m = np.median(x)
    s = np.median(np.abs(x - m))
    u = (x - m) / ((c * s) + epsilon)
    i = np.abs(u) > 1
    w = (1 - u ** 2) ** 2
    w[i] = 0
    tbi = np.sum(w * x) / np.sum(w)
    return tbi


def roku(x):
    # roku
    tbi = tukey_biweight(x)
    vector_p = np.abs(x - tbi)
    h = entropy(vector_p)
    return h


def rankdata(a, method='average'):
    # function for specs
    arr = np.ravel(np.asarray(a))
    algo = 'mergesort' if method == 'ordinal' else 'quicksort'
    sorter = np.argsort(arr, kind=algo)
    inv = np.empty(sorter.size, dtype=np.intp)
    inv[sorter] = np.arange(sorter.size, dtype=np.intp)
    arr = arr[sorter]
    obs = np.r_[True, arr[1:] != arr[:-1]]
    dense = obs.cumsum()[inv]
    # cumulative counts of each unique value
    count = np.r_[np.nonzero(obs)[0], len(obs)]
    # average method
    return .5 * (count[dense] + count[dense - 1] + 1)


def mannwhitneyu(x, y):
    # function for specs
    # mann whitney u test
    x = np.asarray(x)
    y = np.asarray(y)
    n1 = len(x)
    n2 = len(y)
    ranked = rankdata(np.concatenate((x, y)))
    rankx = ranked[0:n1]  # get the x-ranks
    u1 = n1 * n2 + (n1 * (n1 + 1)) / 2.0 - np.sum(rankx, axis=0)  # calc U for x
    u2 = n1 * n2 - u1  # remainder is U for y
    return (u2)


def variance_pi(x1, x2):
    # function for specs
    n1 = len(x1)
    n2 = len(x2)
    wc_stat = mannwhitneyu(x1, x2)
    pi = wc_stat / (n1 * n2)
    return (pi)
