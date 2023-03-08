import numpy as np
import pandas as pd
from scores.intermediate import js_distance, variance_pi

# tissue specific functions

def TSI_spec(x, **kwargs):
    # tissue specificity index
    if not np.any(x):
        return np.array([0.0])
    else:
        n = len(x)
        if n == 1:
            return np.array([0.0])
        else:
            TSI = x / np.sum(x)
            return TSI


def zscore(x, **kwargs):
    # Z-Score
    transform = kwargs.pop('transform', True)
    n = len(x)
    if n == 1:
        return np.array([0.0])
    else:
        std = np.std(x, ddof=1)
        if std == 0:
            return np.zeros(n)
        else:
            zs = (x - np.mean(x)) / std
            if transform:
                max_zs = (n - 1) / np.sqrt(n)
                zs_transformed = (zs + max_zs) / (2 * max_zs)  # normalize to [0,1] interval
                return zs_transformed
            else:
                return zs


def SPM(x, **kwargs):
    # Specificity Measure
    n = len(x)
    if n == 1:
        return np.array([0.0])
    else:
        spm_vector = []
        for i in range(n):
            if x[i] == 0:
                spm_vector.append(0)
            else:
                spm_i = (x[i] ** 2) / (np.linalg.norm(x) * x[i])
                spm_vector.append(spm_i)
        return np.array(spm_vector)


def JSS(x, **kwargs):
    # Jensen-Shannon Specificity
    n = len(x)
    if n == 1:
        return np.array([0.0])
    else:
        js_vector = []
        for i in range(n):
            if x[i] == 0:
                js_vector.append(0)
            else:
                vector_i = np.zeros(n)
                vector_i[i] = x[i]
                js_i = 1 - js_distance(x, vector_i)
                js_vector.append(js_i)
        return np.array(js_vector)


def SPECS(tissues, genes):
    # SPECS

    samples = tuple(tissues.readline().rstrip().split('\t'))
    samples = samples[1:]

    rnas = genes.read().split('\n')
    rnas = list(filter(None, rnas[1:]))

    diseases = tuple([sample.split('_', 1)[0] for sample in samples])
    diseases = np.array(diseases)

    table = pd.Series(diseases)
    counts = table.value_counts(sort=False)
    tot_dis = len(table)
    dis = len(counts)

    unique_dis = []
    seen_dis = set()
    for i in diseases:
        if i not in seen_dis:
            unique_dis.append(i)
            seen_dis.add(i)
    print(unique_dis)

    # calculate weight factor, in this case equally weighted
    pk = 1 / (len(unique_dis) - 1)

    mw_RNA = np.zeros((len(rnas), dis))
    s_RNA = np.zeros((len(rnas), dis))
    transcript = tissues.readline()
    i = 0

    # calculation of score for each transcript
    while transcript:
        start = time.time()
        print(i)
        tmp_RNA = transcript.split('\t')
        tmp_RNA = tmp_RNA[1:]
        tmp_RNA = np.array([float(value) if value != '' else 0.0 for value in tmp_RNA])
        print(tmp_RNA)
        j = 0
        # run through next loop for each sample d
        for sel_dis in unique_dis:
            # print(sel_dis)
            tmp_RNA_selected = tmp_RNA[diseases == sel_dis]
            # print(tmp_RNA_selected)
            unique_dis_non = [unique_dis[k] for k in range(len(unique_dis)) if unique_dis[k] != sel_dis]
            ptot_dis = 0
            var_list = []
            p_list = []
            # compare with each sample not d (k)
            cov_matrix = np.zeros(shape=(len(unique_dis_non), len(unique_dis_non)))
            k = 0
            for comp_dis in unique_dis_non:
                # print(comp_dis)
                l = 0
                tmp_RNA_comp = tmp_RNA[diseases == comp_dis]
                wc_stat_kd = variance_pi(tmp_RNA_selected, tmp_RNA_comp)
                # (wc_stat_kd,var_pi_kd) = variance_pi(tmp_RNA_selected,tmp_RNA_comp)
                # print(var_pi_kd)
                ptot_dis = ptot_dis + wc_stat_kd * pk
                # cov_matrix[k,k] = var_pi_kd
                p_list.append(wc_stat_kd)
            mw_RNA[i, j] = ptot_dis
            j += 1
        transcript = tissues.readline()
        i += 1
        end = time.time()
        print ("Time elapsed:", end - start)

    # writing output
    SPECS_score = pd.DataFrame(mw_RNA, index=rnas, columns=unique_dis)

    return SPECS_score

