from __future__ import print_function, division
from collections import Counter

from keras.models import Sequential, Model
from keras.layers import Dense, Dropout, BatchNormalization, concatenate, Input, LSTM, GRU, Conv1D, add, Flatten
from keras.layers.advanced_activations import PReLU
from keras.optimizers import Nadam
from keras.callbacks import ReduceLROnPlateau

from .preprocess import *
from .eval import *


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


def rnn_model(shape, output, h_units, rnn_type):
    # h_units[0] - number of units in RNN
    # rnn_type = ["lstm", "gru", "bilstm", "bigru"]
    model = Sequential()
    
    if rnn_type.find("gru") != -1:  # TODO: change this to something more sane
        rnn_layer =  GRU(h_units[0], kernel_initializer="he_normal", recurrent_initializer="he_normal", 
                       implementation=2, bias_initializer="he_normal", dropout=.2, recurrent_dropout=.2,
                       unroll=True, input_shape=shape)
    elif rnn_type.find("lstm") != -1:
        rnn_layer = LSTM(h_units[0], kernel_initializer="he_normal", recurrent_initializer="he_normal", 
                       implementation=2, bias_initializer="he_normal", dropout=.2, recurrent_dropout=.2,
                       unroll=True, input_shape=shape)
    else:
        print("Can't find neither GRU not LSTM")
        return 0
    
    model.add(rnn_layer)
    model.add(BatchNormalization())
    model.add(PReLU())
    
    for num in h_units[1:]:
        model.add(Dense(num))
        model.add(BatchNormalization())
        model.add(PReLU())
        model.add(Dropout(.3))
        
    model.add(Dense(output))
    model.add(PReLU())
    
    # model.compile(optimizer="nadam", loss=make_mean_sample_error(max_pos))
    model.compile(optimizer="nadam", loss="mse")
    
    return model


def cnn_pos_model(shape, output, h_units):
    def res_block(prev_layer, shape):
        branch = BatchNormalization()(prev_layer)
        branch = PReLU()(branch)
        branch = Conv1D(h_units[1], 1, kernel_initializer="he_normal")(branch)
        
        branch = BatchNormalization()(branch)
        branch = PReLU()(branch)
        branch = Conv1D(shape[1], 1, kernel_initializer="he_normal")(branch)
        
        return add([prev_layer, branch])
    
    # merged = []
    # merged.append(concatenate([pep_in, len_in]))
    
    # pep_br = Dense(h_units[0])(merged[-1])
    # pep_br = BatchNormalization()(pep_br)
    # pep_br = PReLU()(pep_br)
    # pep_br = Dropout(.3)(pep_br)
    
    # for num in h_units[1:]:
    #     merged.append(concatenate([pep_br, len_in]))
    #     pep_br = Dense(num)(merged[-1])
    #     pep_br = BatchNormalization()(pep_br)
    #     pep_br = PReLU()(pep_br)
    #     pep_br = Dropout(.3)(pep_br)
    
    pep_in = Input(shape)
    pos_in = Input((1,))
    
    pep_branch = res_block(pep_in, shape)
    for ind in range(1, h_units[0]):
        pep_branch = res_block(pep_branch, shape)
        
    pep_branch = Flatten()(pep_branch)
    
    merged = concatenate([pep_branch, pos_in])
    
    merged = Dense(64)(merged)
    merged = BatchNormalization()(merged)
    merged = PReLU()(merged)
    merged = Dropout(.3)(merged)
    
    pred = Dense(output)(merged)
    pred = PReLU()(pred)
    
    model = Model(inputs=[pep_in, pos_in], outputs=pred)
    
    model.compile(optimizer="nadam", loss="mse")
    
    return model


def delta_model(shape, output, h_units):
    def _block(prev_layer, shape):
        branch = BatchNormalization()(prev_layer)
        branch = PReLU()(branch)
        branch = Conv1D(192, 1, kernel_initializer="he_normal")(branch)
        
        branch = BatchNormalization()(branch)
        branch = PReLU()(branch)
        branch = Conv1D(shape[1], 1, kernel_initializer="he_normal")(branch)
        
        return add([prev_layer, branch])
    
    inp_forw = Input(shape = shape)
    inp_back = Input(shape = shape)
    inp_pos = Input(shape = (1,))
    # inp_len = Input((1,)) # one hot encoding for length
    
    shared_model = Sequential()
    shared_model.add(GRU(h_units[0][0], kernel_initializer="he_normal", recurrent_initializer="he_normal", 
                       implementation=2, bias_initializer="he_normal", dropout=.2, recurrent_dropout=.2, 
                       unroll=True, input_shape = shape))
    
    for num in h_units[0][1:]:
        shared_model.add(Dense(num))
        shared_model.add(BatchNormalization())
        shared_model.add(PReLU())
        shared_model.add(Dropout(.3))
    
    diff_forw = shared_model(inp_forw)
    diff_forw = Dense(1)(diff_forw)
    pred_forw = PReLU(name="pred_forw")(diff_forw)
    
    diff_back = shared_model(inp_back)
    diff_back = Dense(1)(diff_back)
    pred_back = PReLU(name="pred_back")(diff_back)
    
    merged = concatenate([pred_forw, pred_back])
    
    for num in h_units[1]:
        merged = concatenate([merged, inp_pos])
        merged = Dense(num)(merged)
        merged = BatchNormalization()(merged)
        merged = PReLU()(merged)
        merged = Dropout(.3)(merged)
    
    merged = Dense(1)(merged)
    pred_coord = PReLU(name="pred_final")(merged)
    
    model = Model(inputs=[inp_forw, inp_back, inp_pos], outputs=[pred_forw, pred_back, pred_coord])
    
    model.compile(optimizer="nadam", loss="mse")
    
    return model


def train_model(max_pos, n_clust, coord, layers, 
                left_window, right_window, 
                n_epochs, 
                hist, model_list, df_loss,
                how = "no", scale = "no", 
                features = "onehot", model_type = "dense", fading = False):
    model_name = model_type + ".l" + str(left_window) + "_r" + str(right_window) + "." + "-".join(map(str, layers)) + "." + how + "_" + scale + "." + features

    if n_clust > 0:
        if fading:
            model_name += ".clust.fade_" + str(n_clust)
        else:
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
        input_shape = (0,)
        if features == "onehot":
            coord_fun = onehot
        elif features == "omega":
            coord_fun = onehot_omega
        elif features == "twohot":
            coord_fun = twohot_omega
        else:
            print("Unknown parameter", features)
            return 0
        
        for_rnn = False
        if model_type in ["gru", "lstm", "cnn_pos"]:
            for_rnn = True
        
        X_can, y_can = coord_fun(df_can, left_window, right_window, max_pos, for_rnn)
        X_cdr, y_cdr = coord_fun(df_cdr, left_window, right_window, max_pos, for_rnn)

        #
        # Prepare to build the model
        #
        if model_type == "dense":
            model_fun = dense_model
            input_shape = (len(CHARS)*(right_window+left_window+1),)
            
        elif model_type == "dense_pos":
            model_fun = dense_pos_model
            input_shape = (len(CHARS)*(right_window+left_window+1),)
            
            # add positions
            X_can = [X_can, np.array([float((x % max_pos) + 1) / max_pos for x in range(X_can.shape[0])])]
            X_cdr = [X_cdr, np.array([float((x % max_pos) + 1) / max_pos for x in range(X_cdr.shape[0])])]
            
        elif model_type == "dense_poslen":
            model_fun = dense_poslen_model
            input_shape = (len(CHARS)*(right_window+left_window+1),)
            
            # add positions and lengths
            X_can = [X_can, np.array([float((x % max_pos) + 1) / max_pos for x in range(X_can.shape[0])]), np.full((X_can.shape[0],1), max_pos)]
            X_cdr = [X_cdr, np.array([float((x % max_pos) + 1) / max_pos for x in range(X_cdr.shape[0])]), np.full((X_cdr.shape[0],1), max_pos)]

        elif model_type in ["gru", "lstm"]:
            model_fun = rnn_model
            input_shape = (right_window+left_window+1, len(CHARS))
        elif model_type in ["cnn_pos"]:
            X_can = [X_can, np.array([float((x % max_pos) + 1) / max_pos for x in range(X_can.shape[0])])]
            X_cdr = [X_cdr, np.array([float((x % max_pos) + 1) / max_pos for x in range(X_cdr.shape[0])])]
            
            model_fun = cnn_pos_model
            input_shape = (right_window+left_window+1, len(CHARS))
        else:
            print("Unknown parameter", coord_fun)
            return 0
        
        if features == "twohot":
            input_shape = (right_window+left_window, 2*len(CHARS))
        
        #
        # Build the model
        #
        if model_type in ["dense", "dense_pos", "dense_poslen"]:
            model = model_fun(input_shape, 1, layers)
        elif model_type in ["gru", "lstm", "cnn_pos"]:
            model = model_fun(input_shape, 1, layers)

        reduce_lr = ReduceLROnPlateau(monitor='val_loss', factor=0.2, patience=3, cooldown=1, min_lr=0.001)
        
        if n_clust == 0:
            hist_obj = model.fit(X_can, y_can, batch_size=64, epochs=n_epochs, verbose=0, validation_data=(X_cdr, y_cdr), callbacks=[reduce_lr])
        else:
            labels, labels_cnt, min_cluster, min_cluster_size = [], -1, -1, -1
            if model_type == "dense":
                labels, labels_cnt, min_cluster, min_cluster_size = cluster(y_can.reshape((-1, max_pos)), n_clust)
            labels_new = []
            for x in labels:
                labels_new.extend([x] * max_pos)
            labels = labels_new
            
            if fading:
                weight_vec = np.array([np.log(min_cluster_size) / np.log(labels_cnt[x]) for x in labels])
                weight_vec = np.exp(np.log(weight_vec) / (200 ** .5))

                hist_obj = model.fit(X_can, y_can, sample_weight=weight_vec, batch_size=64, epochs=n_epochs, verbose=0, validation_data=(X_cdr, y_cdr), callbacks=[reduce_lr])
            else:
                weight_vec = np.array([np.log(min_cluster_size) / np.log(labels_cnt[x]) for x in labels])

                hist_obj = model.fit(X_can, y_can, sample_weight=weight_vec, batch_size=64, epochs=n_epochs, verbose=0, validation_data=(X_cdr, y_cdr), callbacks=[reduce_lr])
        
        hist[model_name] = hist_obj
        model_list[model_name] = model
        
        print(hist_obj.history["val_loss"][-1], end="\t")
        
        boot_loss_vec = bootstrap_cdr(model_list[model_name], X_cdr, y_cdr, max_pos)
        df_new = pd.DataFrame({"val_loss": boot_loss_vec, "model": model_name})
        
        print("(", np.mean(boot_loss_vec), ")")
        
        return model, pd.concat([df_loss, df_new])
    else:
        return model_list[model_name], None