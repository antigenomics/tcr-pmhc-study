import keras
from keras.models import Sequential
from keras.layers import Dense, Dropout, BatchNormalization
from keras.layers.advanced_activations import PReLU


# best model out of Dense models
def dense2_win3():
    model = Sequential()
    
    model.add(Dense(64, input_shape=(140,)))
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
    
    return model