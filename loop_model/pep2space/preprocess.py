import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler, MaxAbsScaler, minmax_scale, maxabs_scale


#
# One-hot
#
CHARS = ["A", "L", "R", 'K', 'N', 'M', 'D', 'F', 'C', 'P', 'Q', 'S', 'E', 'T', 'G', 'W', 'H', 'Y', 'I', 'V']
CHAR_INDICES = dict((c, i) for i, c in enumerate(CHARS))
INDICES_CHAR = dict((i, c) for i, c in enumerate(CHARS))

CHAR_ONE_HOT = dict((c, np.zeros((len(CHARS),), dtype=bool)) for c in CHARS)
for c in CHAR_ONE_HOT:
    CHAR_ONE_HOT[c][CHAR_INDICES[c]] = 1
    
    
def window_to_seq(X, left_window, right_window, max_pos):
    return X.reshape((X.shape[0], (left_window+right_window+1), len(CHARS)))


def onehot(df, left_window, right_window, max_pos, for_rnn = False):
    X = np.zeros((len(df)*max_pos, (left_window+right_window+1) * len(CHARS)), dtype=bool)
    y = np.zeros((len(df)*max_pos, 1), dtype=np.float32)
    for seq_i, seq in enumerate(df["sequence"]):
        seq = "X"*left_window + seq + "X"*right_window
        for index, target_pos in enumerate(range(left_window + 1, len(seq) - right_window)):
            target_aa = seq[target_pos]
            for amb_pos, amb_aa in enumerate(seq[target_pos-left_window : target_pos+right_window+1]):
                if amb_aa != "X":
                    X[seq_i*max_pos + index, amb_pos*len(CHARS):(amb_pos+1)*len(CHARS)] = CHAR_ONE_HOT[amb_aa]
            y[seq_i*max_pos + index] = df[[4 + index]].iloc[seq_i]
            
    if for_rnn:
        X = window_to_seq(X, left_window, right_window, max_pos)
            
    return X, y


def onehot_omega(df, left_window, right_window, max_pos, for_rnn = False, add_pos = False, add_len = False):
    X = np.zeros((len(df)*max_pos, (left_window+right_window+1+add_pos+add_len) * len(CHARS)), dtype=bool)
    y = np.zeros((len(df)*max_pos, 1), dtype=np.float32)
    for seq_i, seq in enumerate(df["sequence"]):
        seq = seq[len(seq) - left_window : len(seq)] + seq + seq[:right_window]
        for index, target_pos in enumerate(range(left_window + 1, len(seq) - right_window)):
            target_aa = seq[target_pos]
            for amb_pos, amb_aa in enumerate(seq[target_pos-left_window : target_pos+right_window+1]):
                if amb_aa != "X":
                    X[seq_i*max_pos + index, amb_pos*len(CHARS):(amb_pos+1)*len(CHARS)] = CHAR_ONE_HOT[amb_aa]
                    if add_pos:
                        X[seq_i*max_pos + index, -1] = float(amb_pos) / max_pos
                    # if add_len:
                    #     X[seq_i*max_pos + index, -1] = max_pos # categotical
            y[seq_i*max_pos + index] = df[[4 + index]].iloc[seq_i]
    
    if for_rnn:
        X = window_to_seq(X, left_window, right_window, max_pos)
        
    return X, y


def twohot_omega(df, left_window, right_window, max_pos, for_rnn = True):
    X = np.zeros((len(df)*max_pos, left_window+right_window, len(CHARS)*2), dtype=bool)
    y = np.zeros((len(df)*max_pos, 1), dtype=np.float32)
    for seq_i, seq in enumerate(df["sequence"]):
        seq = seq[len(seq) - left_window : len(seq)] + seq + seq[:right_window]
        for index, target_pos in enumerate(range(left_window + 1, len(seq) - right_window)):
            target_aa = seq[target_pos]
            
            k = 0
            for amb_pos, amb_aa in enumerate(seq[target_pos-left_window : target_pos+right_window+1]):
                if amb_pos + left_window != target_pos:
                    X[seq_i*max_pos + index, amb_pos - k, :] = np.vstack([CHAR_ONE_HOT[amb_aa].reshape(20,1), CHAR_ONE_HOT[target_aa].reshape(20,1)]).reshape((40,))
                else:
                    k += 1
            
            y[seq_i*max_pos + index] = df[[4 + index]].iloc[seq_i]
    return X, y


def scale_data(df, how, scale):
    if how == "col":
        for i in range(4, 16):
            if scale == "mm":
                df.iloc[:,i] = minmax_scale(df.iloc[:,i])
            elif scale == "abs":
                df.iloc[:,i] = maxabs_scale(df.iloc[:,i])
            else:
                print("Unknown parameter", scale)
    elif how == "all":
        if scale == "mm":
            df.iloc[:,range(4, 16)] = minmax_scale(df.iloc[:,range(4, 16)])
        elif scale == "abs":
            df.iloc[:,range(4, 16)] = maxabs_scale(df.iloc[:,range(4, 16)])
        else:
            print("Unknown parameter", scale)
    else:
        print("Unknown parameter", how)
        
        
def abs_to_diff(y, max_pos, rev = False):
    # We suppose that the first position was zero
    ynew = np.zeros(y.shape)
    for i in range(0, y.shape[0], max_pos):
        if not rev:
            ynew[i] = y[i]
            ynew[i+1 : i+max_pos] = y[i+1 : i+max_pos] - y[i : i+max_pos-1]
        else:
            ynew[i : i+max_pos-1] = y[i : i+max_pos-1] - y[i+1 : i+max_pos]
            ynew[i+max_pos-1] = y[i+max_pos-1]
    return ynew
        

def diff_to_abs(y, max_pos, rev = False):
    ynew = np.zeros(y.shape)
    for i in range(0, y.shape[0], max_pos):
        if not rev:
            ynew[i] = y[i]
            for j in range(1, max_pos):
                ynew[i+j] = y[i+j] + ynew[i+j-1]
        else:
            ynew[i+max_pos-1] = y[i+max_pos-1]
            for j in range(max_pos-1, 0, -1):
                ynew[i+j-1] = y[i+j-1] + ynew[i+j]
    return ynew