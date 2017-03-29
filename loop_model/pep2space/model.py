from collections import Counter

from keras.models import Sequential
from keras.layers import Dense, Dropout, BatchNormalization
from keras.layers.advanced_activations import PReLU
from keras.optimizers import Nadam

from sklearn.cluster import MiniBatchKMeans

from preprocess import *


def add_dense(model, h_units, dropout = .3):
    model.add(Dense(h_units))
    model.add(BatchNormalization())
    model.add(PReLU())
    model.add(Dropout(dropout))

              
def dense_model(shape, output, h_units = [256, 128, 64]):
    model = Sequential()
    
    model.add(Dense(h_units[0], input_shape=shape))
    model.add(BatchNormalization())
    model.add(PReLU())
    model.add(Dropout(.3))
    
    for num in h_units[1:]:
        add_dense(model, num)
        
    model.add(Dense(output))
    model.add(PReLU())
    
    model.compile(optimizer="nadam", loss="mse")
    
    return model


def train_model(max_pos, n_clust, coord, layers, 
                left_window, right_window, 
                n_epochs, 
                hist, model_list, 
                how = "no", scale = "no", 
                coord_fun = "onehot"):
    model_name = "left" + str(left_window) + "_right" + str(right_window) + "." + "-".join(map(str, layers)) + "." + how + "_" + scale + "." + coord_fun

    if n_clust > 0:
        model_name += ".clust_" + str(n_clust)
    
    if model_name not in hist:
        print(model_name)
        df_cdr = pd.read_csv("data/cdr_coord_" + coord + ".csv.gz")
        df_can = pd.read_csv("data/can_coord_" + coord + ".csv.gz")
        
        if how in ["col", "abs"]:
            scale_data(df_cdr, how, scale)
            scale_data(df_can, how, scale)

        if coord_fun == "onehot":
            coord_fun = onehot
        elif coord_fun == "omega":
            coord_fun = onehot_omega
        else:
            print("Unknown parameter", coord_fun)
            return 0
        
        X_can, y_can = coord_fun(df_can, left_window, right_window, max_pos)
        X_cdr, y_cdr = coord_fun(df_cdr, left_window, right_window, max_pos)        

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