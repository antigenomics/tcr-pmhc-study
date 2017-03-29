import pandas as pd
import numpy as np
import keras
from keras.models import Sequential
from keras.layers import Dense, Dropout, BatchNormalization
from keras.layers.advanced_activations import PReLU
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import joblib
from keras import optimizers
from keras.callbacks import ReduceLROnPlateau



################################
################################

coord = "x"
df_cdr = pd.read_csv("data/cdr_coord_" + coord + ".csv.gz")
df_put = pd.read_csv("data/put_coord_" + coord + ".csv.gz")
df_can = pd.read_csv("data/can_coord_" + coord + ".csv.gz")

print(df_cdr.shape)
print(df_put.shape)
print(df_can.shape)


################################
################################

MAX_POS=12

chars = ["A", "L", "R", 'K', 'N', 'M', 'D', 'F', 'C', 'P', 'Q', 'S', 'E', 'T', 'G', 'W', 'H', 'Y', 'I', 'V']
char_indices = dict((c, i) for i, c in enumerate(chars))
indices_char = dict((i, c) for i, c in enumerate(chars))

one_hot = dict((c, np.zeros((len(chars),), dtype=bool)) for c in chars)
for c in one_hot:
    one_hot[c][char_indices[c]] = 1
    
    
################################
################################

def to_vec_onehot(df, left_window, right_window):   
    X = np.zeros((len(df)*MAX_POS, (left_window+right_window+1) * len(chars)+1), dtype=bool)
    y = np.zeros((len(df)*MAX_POS, 1), dtype=np.float32)
    for seq_i, seq in enumerate(df["sequence"]):
        seq = "X"*left_window + seq + "X"*right_window
        for index, target_pos in enumerate(range(left_window + 1, len(seq) - right_window)):
            target_aa = seq[target_pos]
            for amb_pos, amb_aa in enumerate(seq[target_pos-left_window : target_pos+right_window+1]):
                if amb_aa != "X":
                    X[seq_i*MAX_POS + index, amb_pos*len(chars):(amb_pos+1)*len(chars)] = one_hot[amb_aa]
            X[seq_i*MAX_POS + index, -1] = target_pos / len(seq)
            y[seq_i*MAX_POS + index] = df[[4 + index]].iloc[seq_i]
    return X, y


def to_vec_kidera(df, left_window, right_window):
    X = np.zeros((len(df)*MAX_POS, (left_window+right_window+1) * 10), dtype=bool)
    y = np.zeros((len(df)*MAX_POS, 1), dtype=np.float32)
    for seq_i, seq in enumerate(df["sequence"]):
        seq = "X"*left_window + seq + "X"*right_window
        for index, target_pos in enumerate(range(left_window + 1, len(seq) - right_window)):
            target_aa = seq[target_pos]
            for amb_pos, amb_aa in enumerate(seq[target_pos-left_window : target_pos+right_window+1]):
                if amb_aa != "X":
                    X[seq_i*MAX_POS + index, amb_pos*10:(amb_pos+1)*10] = kidera[amb_aa]
            y[seq_i*MAX_POS + index] = df[[4 + index]].iloc[seq_i]
    return X, y


################################
################################

X_can, y_can = to_vec_onehot(df_can, 3, 3)
print(X_can.shape)
# X_put, y_put = to_vec_onehot(df_put, 3, 3)
X_put = np.load("X_put.npy")
y_put = np.load("y_put.npy")
print(X_put.shape)
X_cdr, y_cdr = to_vec_onehot(df_cdr, 3, 3)
print(X_cdr.shape)


################################
################################

def add_dense(model, h_units):
    model.add(Dense(h_units))
    model.add(BatchNormalization())
    model.add(PReLU())
    model.add(Dropout(.3))

              
def dense_model(shape, output, h_units = [256, 128, 64]):
    model = Sequential()
    
    model.add(Dense(h_units[0], input_shape=shape))
    model.add(BatchNormalization())
    model.add(PReLU())
    model.add(Dropout(.3))
    
    for num in h_units[1:]: add_dense(model, num)
        
    model.add(Dense(output))
    model.add(PReLU())
    
    model.compile(optimizer="nadam", loss="mse")
    
    return model


from collections import Counter


def to_vec_onehot(df, left_window, right_window):
    X = np.zeros((len(df)*MAX_POS, (left_window+right_window+1) * len(chars)), dtype=bool)
    y = np.zeros((len(df)*MAX_POS, 1), dtype=np.float32)
    for seq_i, seq in enumerate(df["sequence"]):
        seq = "X"*left_window + seq + "X"*right_window
        for index, target_pos in enumerate(range(left_window + 1, len(seq) - right_window)):
            target_aa = seq[target_pos]
            for amb_pos, amb_aa in enumerate(seq[target_pos-left_window : target_pos+right_window+1]):
                if amb_aa != "X":
                    X[seq_i*MAX_POS + index, amb_pos*len(chars):(amb_pos+1)*len(chars)] = one_hot[amb_aa]
            y[seq_i*MAX_POS + index] = df[[4 + index]].iloc[seq_i]
    return X, y


def preprocess(df, how, scale):
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


def train_models(n_clust, coord, layers, left_window, right_window, n_epochs, hist, model_list, how = "no", scale = "no"):
    model_name = "left" + str(left_window) + "_right" + str(right_window) + "." + "-".join(map(str, layers)) + "." + how + "_" + scale

    if n_clust > 0:
        model_name += ".clust_" + str(n_clust)
    
    if model_name not in hist:
        print(model_name)
        df_cdr = pd.read_csv("data/cdr_coord_" + coord + ".csv.gz")
        df_can = pd.read_csv("data/can_coord_" + coord + ".csv.gz")
        
        if how in ["col", "abs"]:
            preprocess(df_cdr, how, scale)
            preprocess(df_can, how, scale)

        X_can, y_can = to_vec_onehot(df_can, left_window, right_window)
        X_cdr, y_cdr = to_vec_onehot(df_cdr, left_window, right_window)        

        model = dense_model((20*(right_window+left_window+1),), 1, layers)

        if n_clust == 0:
            hist_obj = model.fit(X_can, y_can, batch_size=64, epochs=n_epochs, verbose=0, validation_data=(X_cdr, y_cdr))
        else:
            kmeans = MiniBatchKMeans(n_clust)
            kmeans.fit(X_can)
            
            labels = kmeans.predict(X_can)
            labels_cnt = Counter(labels)
            min_cluster, min_cluster_size = min(labels_cnt.items(), key = lambda x: x[1])
            
            weight_vec = np.array([np.log(min_cluster_size) / np.log(labels_cnt[x]) for x in labels])
            
            hist_obj = model.fit(X_can, y_can, sample_weight=weight_vec, batch_size=64, epochs=n_epochs, verbose=0, validation_data=(X_cdr, y_cdr))
        
        hist[model_name] = hist_obj
        model_list[model_name] = model
        
        return model
    else:
        print(model_name, "- return the old model")
        return model_list[model_name]


################################
################################

hist = {}
for left_window in range(8):
    for right_window in range(8):
        train_models("y", left_window, right_window, 2000, hist)


################################
################################

fig, ax = plt.subplots(nrows=1, sharex=True, ncols=2)
fig.set_figwidth(16)
fig.set_figheight(10)

best_models = sorted([(h, np.mean(hist[h].history["val_loss"][-5:])) for h in hist], key=lambda x: x[1])[:8]

for i, (h, _) in enumerate(sorted(best_models)):
    ax[0].plot(np.log2(hist[h].history["loss"][100:]), label=h)
    ax[1].plot(np.log2(hist[h].history["val_loss"][100:]), label=h)


ax[0].set_title("loss")
ax[1].set_title("val")
ax[0].legend()
ax[1].legend()

plt.savefig("loss_y_dense_onehot_8best_2000it.png")