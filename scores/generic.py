import numpy as np
from scipy.stats import moment
from scores.intermediate import entropy, dpm, roku
from scores.tissuespecific import JSS, SPM


def counts(x, **kwargs):
    threshold = kwargs.pop('threshold', 0)
    n = len(x)
    if n <= 1:
        return 0.0
    else:
        cts = np.sum(x > threshold)
        if cts == 0:
            return 0.0
        else:
            cts_transformed = (1 - (cts / n)) * (n / (n - 1))
            return cts_transformed


def tau(x, **kwargs):
    if not np.any(x):
        return 0.0
    else:
        n = len(x)
        if n == 1:
            return 0.0
        else:
            x_hat = x / np.max(x)
            tau = (1 / (n - 1)) * np.sum(1 - x_hat)
            return tau


def gini(x, **kwargs):
    transform = kwargs.pop('transform', True)
    if not np.any(x):
        return 0.0
    else:
        x = np.sort(x)
        n = len(x)
        if n == 1:
            return 0.0
        else:
            index = np.arange(1, n + 1)
            gini_coeff = (np.sum((2 * index - (n + 1)) * x)) / (n * np.sum(x))
            if transform:
                transformed_gini_coefficient = gini_coeff * (n / (n - 1))
                return transformed_gini_coefficient
            else:
                return gini_coeff


def HS(x, **kwargs):
    # shannon entropy
    transform = kwargs.pop('transform', True)
    if not np.any(x):
        return 0.0
    else:
        n = len(x)
        if n == 1:
            return 0.0
        else:
            HS = np.log2(n) - entropy(x)
            if transform:
                HS_transformed = HS / np.log2(n)
                return HS_transformed
            else:
                return HS


def simpson(x, **kwargs):
    # simpson index
    transform = kwargs.pop('transform', True)
    if not np.any(x):
        return 0.0
    else:
        n = len(x)
        if n == 1:
            return 0.0
        else:
            p = x / np.sum(x)
            simpson_index = np.sum(p ** 2)
            if transform:
                min_simpson = 1 / len(x)
                transformed_simpson_index = (simpson_index - min_simpson) / (1 - min_simpson)
                return transformed_simpson_index
            else:
                return simpson_index


def ROKU_spec(x, **kwargs):
    transform = kwargs.pop('transform', True)
    if not np.any(x):
        return 0.0
    else:
        n = len(x)
        if n == 1:
            return 0.0
        else:
            rs = np.log2(n) - roku(x)
            if transform:
                rs_transformed = rs / np.log2(n)
                return rs_transformed
            else:
                return rs


def SPM_DPM(x, **kwargs):
    # specificity measure dispersion
    spm_vector = SPM(x)
    spm_dispersion = dpm(spm_vector)
    return spm_dispersion


def JSS_DPM(x, **kwargs):
    # jensen shannon specificity dispersion
    js_vector = JSS(x)
    js_specificity_dispersion = dpm(js_vector)
    return js_specificity_dispersion


def kurto(x, **kwargs):
    # kurtosis
    n = len(x)
    K = moment(x, moment=4) / (moment(x, moment=2)) ** 2 - 3
    return K
