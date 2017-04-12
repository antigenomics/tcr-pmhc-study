from __future__ import division

import numpy as np
import pandas as pd
import scipy.stats as stats

from sklearn.metrics import mean_squared_error

import keras.backend as K

from .preprocess import *


def get_true_pred_col(filepath, scaler, model, l, r, max_pos, coord_fun = onehot):
    df = pd.read_csv(filepath)
    sc_ls = []
    
    for i in range(4, 4+max_pos):
        sc = scaler()
        df.iloc[:,i] = sc.fit_transform(df.iloc[:,i].values.reshape(-1, 1))
        sc_ls.append(sc)
    
    X, y = coord_fun(df, l, r)

    y_tr = y.reshape((len(df), max_pos))
    for sc_i, i in enumerate(range(4, 4+max_pos)):
        y_tr[:,sc_i] = sc_ls[sc_i].inverse_transform(y_tr[:,sc_i].reshape(-1,1)).reshape((y_tr.shape[0],))
    y_tr = y_tr.reshape((X.shape[0],))
    
    y_pred = model.predict(X)
    y_pred = y_pred.reshape((len(df), max_pos))
    for sc_i, i in enumerate(range(4, 4+max_pos)):
        y_pred[:,sc_i] = sc_ls[sc_i].inverse_transform(y_pred[:,sc_i].reshape(-1,1)).reshape((y_pred.shape[0],))
    y_pred = y_pred.reshape((X.shape[0],))
    
    return y_tr, y_pred


def tr_pred_col(filepath, scaler, model, l, r, coord_fun = onehot):
    y_tr, y_pr = get_true_pred_col(filepath, scaler, model, l, r, max_pos, coord_fun)
    return mean_squared_error(y_tr, y_pr)


def get_true_pred_all(filepath, scaler, model, l, r, max_pos, coord_fun = onehot):
    df = pd.read_csv(filepath)
    df.iloc[:,range(4, 4+max_pos)] = scaler.fit_transform(df.iloc[:,range(4, 4+max_pos)])
    X, y = coord_fun(df, l, r)

    y_tr = y.reshape((len(df), max_pos))
    y_tr = scaler.inverse_transform(y_tr)
    y_tr = y_tr.reshape((X.shape[0],))
    
    y_pred = model.predict(X)
    y_pred = y_pred.reshape((len(df), max_pos))
    y_pred = scaler.inverse_transform(y_pred)
    y_pred = y_pred.reshape((X.shape[0],))
    
    return y_tr, y_pred


def tr_pred_all(filepath, scaler, model, coord_fun = onehot):
    y_tr, y_pr = get_true_pred_all(filepath, scaler, model, l, r, max_pos, coord_fun)
    return mean_squared_error(y_tr, y_pr)


def mean_sample_error(y_true, y_pred, max_pos):
    return ((y_true.reshape((-1, max_pos)).sum(1) - y_pred.reshape((-1, max_pos)).sum(1)) ** 2).mean()


# use decorators later
def make_mean_sample_error(max_pos):
    def mse(y_true, y_pred):
        return K.mean(K.square(K.sum(K.reshape(y_true, (-1, max_pos)), 1) - K.sum(K.reshape(y_pred, (-1, max_pos)), 1)))
        
    return mse


def bootstrap_cdr(model, X, y, max_pos, n = 1000, proc_y = lambda x: x):
    if type(X) is list:
        n_objects = X[0].shape[0] // max_pos
    else:
        n_objects = X.shape[0] // max_pos
    res = []
    for i in range(n):
        inds = np.random.randint(0, n_objects, n_objects)
        to_take = []
        for ind in inds:
            to_take.extend(range(ind * max_pos, (ind + 1) * max_pos))
            
        ynew = y[to_take]
        if type(X) is list:
            Xnew = [x[to_take] for x in X]
        else:
            Xnew = X[to_take]
        res.append(mean_sample_error(model.predict(Xnew), ynew, max_pos))
    return res