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

model_name = "dense2_win3_one"
model = Sequential()
model.add(Dense(128, input_shape=(140,)))
model.add(BatchNormalization())
model.add(Dropout(.3))
model.add(PReLU())

model.add(Dense(64))
model.add(BatchNormalization())
model.add(Dropout(.3))
model.add(PReLU())

# model.add(Dense(64))
# model.add(BatchNormalization())
# model.add(Dropout(.3))
# model.add(PReLU())

model.add(Dense(1))
model.add(PReLU())

# opt=optimizers.Nadam(lr=0.004)
model.compile(optimizer="nadam", loss="mse")


################################
################################

def train_models(coord, left_window, right_window, n_epochs, hist):
    model_name = "left" + str(left_window) + "_right" + str(right_window)
    if model_name not in hist:
        print(model_name)
        df_cdr = pd.read_csv("data/cdr_coord_" + coord + ".csv.gz")
        df_can = pd.read_csv("data/can_coord_" + coord + ".csv.gz")

        X_can, y_can = to_vec_onehot(df_can, left_window, right_window)
        X_cdr, y_cdr = to_vec_onehot(df_cdr, left_window, right_window)

        model = Sequential()
        model.add(Dense(128, input_shape=(20*(right_window+left_window+1)+1,)))
        model.add(BatchNormalization())
        model.add(Dropout(.3))
        model.add(PReLU())

        model.add(Dense(64))
        model.add(BatchNormalization())
        model.add(Dropout(.3))
        model.add(PReLU())

        model.add(Dense(1))
        model.add(PReLU())
        model.compile(optimizer="nadam", loss="mse")

        hist[model_name] = model.fit(X_can, y_can, batch_size=128, epochs=n_epochs, verbose=0, validation_data=(X_cdr, y_cdr))


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