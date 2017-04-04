from __future__ import print_function, division
from collections import Counter

from keras.models import Sequential, Model
from keras.layers import Dense, Dropout, BatchNormalization, concatenate, Input, LSTM, GRU
from keras.layers.advanced_activations import PReLU
from keras.optimizers import Nadam
from keras.callbacks import ReduceLROnPlateau

from sklearn.cluster import MiniBatchKMeans

from .preprocess import *


def dense_model(shape, output, h_units = [256, 128, 64]):  
    model = Sequential()
    
    model.add(Dense(h_units[0], input_shape=shape))
    model.add(BatchNormalization())
    model.add(PReLU())
    model.add(Dropout(.3))
    
    for num in h_units[1:]:
        model.add(Dense(num))
        model.add(BatchNormalization())
        model.add(PReLU())
        model.add(Dropout(.3))
        
    model.add(Dense(output))
    model.add(PReLU())
    
    model.compile(optimizer="nadam", loss="mse")
    
    return model


def dense_pos_model(shape, output, h_units):
    pep_in = Input(shape)
    len_in = Input((1,))
    
    merged = []
    merged.append(concatenate([pep_in, len_in]))
    
    pep_br = Dense(h_units[0])(merged[-1])
    pep_br = BatchNormalization()(pep_br)
    pep_br = PReLU()(pep_br)
    pep_br = Dropout(.3)(pep_br)
    
    for num in h_units[1:]:
        merged.append(concatenate([pep_br, len_in]))
        pep_br = Dense(num)(merged[-1])
        pep_br = BatchNormalization()(pep_br)
        pep_br = PReLU()(pep_br)
        pep_br = Dropout(.3)(pep_br)
    
    merged.append(concatenate([pep_br, len_in]))
    pep_br = Dense(output)(merged[-1])
    pred = PReLU()(pep_br)
    
    model = Model(inputs=[pep_in, len_in], outputs=pred)
    
    model.compile(optimizer="nadam", loss="mse")
    
    return model


def dense_poslen_model(shape, output, h_units):
    pep_in = Input(shape)
    pos_in = Input((1,))
    len_in = Input((1,))
    
    merged = []
    merged.append(concatenate([pep_in, pos_in, len_in]))
    
    pep_br = Dense(h_units[0])(merged[-1])
    pep_br = BatchNormalization()(pep_br)
    pep_br = PReLU()(pep_br)
    pep_br = Dropout(.3)(pep_br)
    
    for num in h_units[1:]:
        merged.append(concatenate([pep_br, pos_in, len_in]))
        pep_br = Dense(num)(merged[-1])
        pep_br = BatchNormalization()(pep_br)
        pep_br = PReLU()(pep_br)
        pep_br = Dropout(.3)(pep_br)
    
    merged.append(concatenate([pep_br, pos_in, len_in]))
    pep_br = Dense(output)(merged[-1])
    pred = PReLU()(pep_br)
    
    model = Model(inputs=[pep_in, pos_in, len_in], outputs=pred)
    
    model.compile(optimizer="nadam", loss="mse")
    
    return model


def rnn_model(shape, output, h_units, rnn_type = "gru"):
    # h_units[0] - number of units in RNN
    # rnn_type = ["lstm", "gru", "bilstm", "bigru"]
    model = Sequential()
    
    if rnn_type.find("gru") != -1:  # TODO: change this to something more sane
        rnn_layer =  GRU(h_units[0], kernel_initializer="he_normal", recurrent_initializer="he_normal", 
                       implementation=2, bias_initializer="he_normal",
                       dropout=.2, recurrent_dropout=.2,
                       unroll=True, input_shape=shape)
    elif rnn_type.find("lstm") != -1:
        rnn_layer = LSTM(h_units[0], kernel_initializer="he_normal", recurrent_initializer="he_normal", 
                       implementation=2, bias_initializer="he_normal",
                       dropout=.2, recurrent_dropout=.2,
                       unroll=True, input_shape=shape)
    else:
        print("Can't find neither GRU not LSTM")
        return 0
    
    model.add(rnn_layer)
    
    for num in h_units[1:]:
        model.add(Dense(num))
        model.add(BatchNormalization())
        model.add(PReLU())
        model.add(Dropout(.3))
        
    model.add(Dense(output))
    model.add(PReLU())
    
    model.compile(optimizer="nadam", loss="mse")
    
    return model


def train_model(max_pos, n_clust, coord, layers, 
                left_window, right_window, 
                n_epochs, 
                hist, model_list, 
                how = "no", scale = "no", 
                features = "onehot", model_type = "dense"):
    model_name = model_type + ".l" + str(left_window) + "_r" + str(right_window) + "." + "-".join(map(str, layers)) + "." + how + "_" + scale + "." + features

    if n_clust > 0:
        model_name += ".clust_" + str(n_clust)
    
    if model_name not in hist:
        print(model_name, end = "\t")
        
        #
        # Load the data
        #
        df_cdr = pd.read_csv("data/cdr_coord_" + coord + ".csv.gz")
        df_can = pd.read_csv("data/can_coord_" + coord + ".csv.gz")
        
        #
        # Scale the data
        #
        if how in ["col", "abs"]:
            scale_data(df_cdr, how, scale)
            scale_data(df_can, how, scale)

        #
        # Extract feactures
        #
        if features == "onehot":
            coord_fun = onehot
        elif features == "omega":
            coord_fun = onehot_omega
        else:
            print("Unknown parameter", features)
            return 0
        
        X_can, y_can = coord_fun(df_can, left_window, right_window, max_pos)
        X_cdr, y_cdr = coord_fun(df_cdr, left_window, right_window, max_pos)        

        #
        # Prepare to build the model
        #
        if model_type == "dense":
            model_fun = dense_model
            
        elif model_type == "dense_pos":
            model_fun = dense_pos_model
            # add positions
            X_can = [X_can, np.array([float((x % max_pos) + 1) / max_pos for x in range(X_can.shape[0])])]
            X_cdr = [X_cdr, np.array([float((x % max_pos) + 1) / max_pos for x in range(X_cdr.shape[0])])]
            
        elif model_type == "dense_poslen":
            model_fun = dense_poslen_model
            # add positions and lengths
            X_can = [X_can, np.array([float((x % max_pos) + 1) / max_pos for x in range(X_can.shape[0])]), np.full((X_can.shape[0],1), max_pos)]
            X_cdr = [X_cdr, np.array([float((x % max_pos) + 1) / max_pos for x in range(X_cdr.shape[0])]), np.full((X_cdr.shape[0],1), max_pos)]

        elif model_type in ["gru", "lstm"]:
            model_fun = rnn_model
            X_can = window_to_seq(X_can, left_window, right_window, max_pos)
            X_cdr = window_to_seq(X_cdr, left_window, right_window, max_pos)
        else:
            print("Unknown parameter", coord_fun)
            return 0
        
        #
        # Build the model
        #
        if model_type in ["dense", "dense_pos", "dense_poslen"]:
            model = model_fun((len(CHARS)*(right_window+left_window+1),), 1, layers)
        elif model_type in ["gru", "lstm"]:
            model = model_fun((right_window+left_window+1, len(CHARS)), 1, layers, model_type)

        reduce_lr = ReduceLROnPlateau(monitor='val_loss', factor=0.2, patience=3, cooldown=1, min_lr=0.0005)
        
        if n_clust == 0:
            hist_obj = model.fit(X_can, y_can, batch_size=64, epochs=n_epochs, verbose=0, validation_data=(X_cdr, y_cdr), callbacks=[reduce_lr])
        else:
            kmeans = MiniBatchKMeans(n_clust, batch_size=1000, n_init=10)
            
            if model_type == "dense":
                kmeans.fit(X_can)
                labels = kmeans.predict(X_can)
            else:
                kmeans.fit(X_can[0])
                labels = kmeans.predict(X_can[0])
                
            labels_cnt = Counter(labels)
            min_cluster, min_cluster_size = min(labels_cnt.items(), key = lambda x: x[1])
            
            weight_vec = np.array([np.log(min_cluster_size) / np.log(labels_cnt[x]) for x in labels])
            
            hist_obj = model.fit(X_can, y_can, sample_weight=weight_vec, batch_size=64, epochs=n_epochs, verbose=0, validation_data=(X_cdr, y_cdr), callbacks=[reduce_lr])
        
        hist[model_name] = hist_obj
        model_list[model_name] = model
        
        print(hist_obj.history["val_loss"][-1])
        
        return model
    else:
        print(model_name, "- return the old model")
        return model_list[model_name]