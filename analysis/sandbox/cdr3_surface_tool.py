
# coding: utf-8

# In[201]:

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import time
from scipy.optimize import minimize

alphabet = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', \
			'G', 'H', 'I', 'L', 'K', 'M', 'F', \
			'P', 'S', 'T', 'W', 'Y', 'V']
			
aa_props = pd.read_csv('amino_acid_properties.txt', sep='\t', index_col=0)

PROP_NUM = 3
hydfob_unnormed = aa_props.hydrophobicity                      
charge_unnormed = aa_props.charge
volume_unnormed = aa_props.volume

def normalize_series(series):
    mx = max(series)
    mn = min(series)
    fun = lambda x: (2.0*x-(mx+mn))/(mx-mn)
    return series.apply(fun)

hydfob = normalize_series(hydfob_unnormed)
charge = normalize_series(charge_unnormed)
volume = normalize_series(volume_unnormed)


# In[2]:

# Fetch and correct peptide table

def discretize(x):
    if x == 'immunogenic':
        return 1
    if x == 'non-immunogenic':
        return 0
    return np.nan

def check_amino_acids(x):
    return all(char in alphabet for char in x)

# In[3]:

# Convert peptides into arrays of properties

def get_surface(peptide):
    surf = []
    for aa in peptide:
        surf.append((hydfob[aa], charge[aa], volume[aa]))
    return tuple(surf)

# In[ ]:

# Functions for calculating score between two peptides
#
# Calculates score for given alignment (shift) of two peptides
# whithout summing the contributions from each amino acid property
def helper(surf1, surf2, shift):
    score = [0] * PROP_NUM
    j = 0
    while j < PROP_NUM:
        i = 0
        while i < len(surf2):
            a1 = surf1[i + shift]
            a2 = surf2[i]
            score[j] += abs(a1[j] - a2[j])
            if i < len(surf2) - 1:
                a11 = surf1[i + shift + 1]
                a21 = surf2[i + 1]
                score[j] += abs(a11[j] - a1[j] - a21[j] + a2[j])
            i += 1
        score[j] /= i
        j += 1
    return score

# Gets minimal score among all alignments of two given peptides
# whithout choosing the minimal one among all alingments of two peptides
def get_surface_distance_expanded(surf1, surf2):
    if (len(surf1) < len(surf2)):
        buf = surf1
        surf1 = surf2
        surf2 = buf
    shift_max = len(surf1) - len(surf2) + 1
    scores = [helper(surf1, surf2, shift) for shift in xrange(shift_max)]
    return scores

# Gets expanded scores for set vs set --> np.array[][]   
def get_set_vs_set(a1, a2):
	res = np.array([[0]*len(a2)]*len(a1), dtype=object)
	for i in xrange(len(a1)):
		for j in xrange(len(a2)):
		    res[i][j] = get_surface_distance_expanded(a1[i], a2[j])
	return res

# In[258]:

def get_score_using_weights(expanded_scores, weights=[1.,  1.]):
    minx = np.inf
    for i in xrange(len(expanded_scores)):
        score = expanded_scores[i][0]
        for j in xrange(len(weights)):
            score += weights[j] * expanded_scores[i][j + 1]
        if score < minx:
            minx = score
    return minx

# In[180]:

def silhouette(self_vs_self, self_vs_other):
    single_point_silhouette_sum = 0
    for i in xrange(self_vs_self.shape[0]):
        self_mean = np.mean(self_vs_self[i])
        other_mean = np.mean(self_vs_other[i])
        single_point_silhouette_sum += (other_mean - self_mean) / max(self_mean, other_mean)
    return single_point_silhouette_sum / self_vs_self.shape[0]


# In[306]:

# just for testing
def minima(*a):
    try:
        g = -600.0 / max(map(abs, a))
    except ZeroDivisionError:
        return a[0]
    den = sum([np.exp(g * x) for x in a])
    num = sum([x * np.exp(g * x) for x in a])
    return num / den

# check 
def check_num(x):
    return (np.isnan(x) or np.isinf(x))

# returns index of min value in the list l
def index_min(l):
    i = 1
    pointer = 0
    minim = l[0]
    while i < len(l):
        if l[i] < minim:
            minim = l[i]
            pointer = i
        i += 1
    return pointer

# returns derivative of min(a[]) function (used when calculating gradient of silhouette)
def min_deriv(a, a_deriv):
    return a_deriv[index_min(a)]

# returns derivative of max(a[]) function
def max_deriv(a, a_deriv):
    i = 1
    pointer = 0
    maxim = a[0]
    while i < len(a):
        if a[i] > maxim:
            maxim = a[i]
            pointer = i
        i += 1
    return a_deriv[pointer]

# tau and psi are lists of average values of distance from each antigen in one cluster to another cluster
# output of this function (psi or tau) depends on the input array 'self_vs_smth'
def tau_psi_calc(weights, self_vs_smth):
    abs_weights = map(np.abs, weights)
    size = len(self_vs_smth[0])
    res  = np.array([0.] * len(self_vs_smth))
    for i in range(len(self_vs_smth)):
        for j in range(size):
            vals = [x[0] + sum([abs_weights[m] * x[m + 1] for m in range(PROP_NUM - 1)]) for x in self_vs_smth[i][j]]
            res[i] += min(vals)
        res[i] /= size
    return res

# returns derivative of tau or psi
def tau_psi_deriv(weights, self_vs_smth):
    abs_weights = map(np.abs, weights)
    sign_weights = map(np.sign, weights)
    size = len(self_vs_smth[0])
    res  = np.array([[0.] * len(self_vs_smth)] * (PROP_NUM - 1), dtype=float)
    for i in range(len(self_vs_smth)):
        for j in range(size):
            vals = [x[0] + sum([abs_weights[m] * x[m + 1] for m in range(PROP_NUM - 1)]) for x in self_vs_smth[i][j]]
            id_min = index_min(vals)
            weight_id = 0
            while weight_id < PROP_NUM - 1:
                min_der = sign_weights[weight_id] * self_vs_smth[i][j][id_min][weight_id + 1]
                res[weight_id, i] = res[weight_id, i] + min_der
                weight_id += 1
        weight_id = 0
        while weight_id < PROP_NUM - 1:
            res[weight_id, i] = res[weight_id, i] / size
            weight_id += 1
    return res
 
# negative gradient of silhouette
def silhouette_grad_negative(weights, self_vs_self, self_vs_other):
    tau = tau_psi_calc(weights, self_vs_other)
    psi = tau_psi_calc(weights, self_vs_self)
    tau_deriv = tau_psi_deriv(weights, self_vs_other)
    psi_deriv = tau_psi_deriv(weights, self_vs_self)
    grad = []
    for weight_id in range(PROP_NUM - 1):
        sil_deriv = i = 0
        tau_deriv_cur = tau_deriv[weight_id]
        psi_deriv_cur = psi_deriv[weight_id]
        while i < len(self_vs_self):
            maxima = max(tau[i], psi[i])
            sil_deriv += (((tau_deriv_cur[i] - psi_deriv_cur[i]) * maxima - (tau[i] - psi[i]) * max_deriv([tau[i], psi[i]], [tau_deriv_cur[i], psi_deriv_cur[i]])) / (maxima * maxima))
            i += 1
        grad.append(-sil_deriv / i)
    return np.array(grad)

# negative value of silhouette (we need maximum of silhouette so we minimize its negative value)
def silhouette_value_negative(weights, self_vs_self, self_vs_other):
    tau = tau_psi_calc(weights, self_vs_other)
    psi = tau_psi_calc(weights, self_vs_self)
    return (-np.mean([(i - j)/max(i, j) for i, j in zip(tau, psi)]))

# In[334]:

# performs gradient descent
def gradient_descent(weights, self_vs_self, self_vs_other, silouette_vals, weights_vals, alpha=1000., max_iter=50, iter_id=0):
    sil_val = silhouette_value_negative(weights, self_vs_self, self_vs_other)
    rel_difference = 1.
    if iter_id:
        rel_difference = abs((sil_val - silouette_vals[-1]) / silouette_vals[-1])
    silouette_vals.append(sil_val)
    weights_vals.append(weights)
    if iter_id < max_iter and (rel_difference > 0.00001 or np.isnan(rel_difference)):
        new_weights = weights - (alpha * silhouette_grad_negative(weights, self_vs_self, self_vs_other))
        gradient_descent(new_weights, self_vs_self, self_vs_other, silouette_vals, weights_vals, alpha, max_iter, iter_id+1)
    else:
        return weights
