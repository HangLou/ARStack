import numpy as np
import pandas as pd
from keras.models import Sequential
from keras.layers import LSTM, Dense, Dropout, Embedding, Activation, Bidirectional
from keras.optimizers import Adam, SGD
from keras import utils
from sklearn.metrics import mean_squared_error as mse
import tensorflow as tf
from keras import backend as B
import os
from collections import Counter
import traceback
import random

numbers = np.arange(0, 1000, 100)
numbers = dict(zip(numbers, numbers))


def get_data(dataset: int):

    csv = pd.read_csv(
        'data/input/processed/MSA_Dataset'+str(dataset)+'.csv')
    #del csv['Unnamed: 0']
    return csv


def transform(y_log):  # no transformations for paper
    #    y_log = 1/(1+np.exp(-y_log))#np.log(abs(y))
    #    y_log=-1*np.log10(y_log)
    return y_log


def transform_back(yp):
    #    yp=10**(-1*yp)
    #    yp= np.log(yp/(1-yp))#np.exp(yp)
    return yp


def get_parameters():
    lr = 0.01  # 0.01 fr paper
    epcs = 3  # 3 for paper
    bsize = 150  # 150 for paper
    return lr, epcs, bsize


def predict(X, y):

    global record
    record = []
    proteins = set(X['pdb'])
    win = 0
    coords = []
    ii = 0
    print(len(proteins))
    for prot in list(proteins):
        ii += 1
        print('number: ', ii)
        try:
            Xts = X.loc[X['pdb'] == prot]
            print(len(Xts))
            Xtr = X.drop(Xts.index)

            ytr = y.loc[Xtr.index]

            Xtr = Xtr.sample(frac=1, random_state=9438)
            ytr = ytr.sample(frac=1, random_state=9438)

            ytr = np.array(ytr)
            yts = np.array(y.loc[Xts.index])

            yts = yts.reshape(-1, 1)
            ytr = ytr.reshape(-1, 1)

            Xtr = Xtr.iloc[:, :-1]
            Xts = Xts.iloc[:, :-1]

            model = []
            model = Sequential()
            model.add(Dense(200, activation='relu',
                      input_shape=(Xtr.shape[1],)))
            model.add(Dropout(0.2))
            model.add(Dense(200, activation='relu'))
            model.add(Dropout(0.2))
            model.add(Dense(100, activation='softmax'))
            # model.add(Dropout(0.2))
            model.add(Dense(1))

            lr, epcs, bsize = get_parameters()

            adam = Adam(learning_rate=0.001)
            sgd = SGD(learning_rate=lr)

            model.compile(loss='mae', optimizer=sgd, metrics=['mse'])

           # model = SVR()  # RandomForestRegressor(n_estimators =100)
            # ,epochs=epcs,verbose=0, batch_size = bsize) # 4, 100 works well. Use at least 6 epochs with the settings above
            #model.fit(Xtr, ytr)

            yp = model.predict(Xts)
            yp = [i[0] for i in yp]

            rmse_pred = (mse(yts, yp)**0.5)
            rmse_alpha = (mse(yts, np.zeros(len(yts)))**0.5)
            rmse_pred = np.round(rmse_pred, 5)
            rmse_alpha = np.round(rmse_alpha, 5)

            print('alpha: ', rmse_alpha)
            print('keras', rmse_pred)
            if rmse_pred < rmse_alpha:
                win += 1
                print('wins: ', win)

            record.append([rmse_alpha, rmse_pred])  # recording RMSEs
            # recording new magnitudes (not coordinates, I will convert these in the next script)
            coords.append([prot]+list(yp))

        except:
            print('error: ', prot)
            print(traceback.print_exc())

    return coords, record


if __name__ == '__main__':
    dataset = 1

    global record
    record = []
    csv = get_data(dataset=1)
    df = csv
    df = df[0::3]
    X = df.iloc[:, 6:]

    X['pdb'] = df['0']
    y = df['2']
    coords, record = predict(X, y)

    df2 = pd.DataFrame(coords)
    df2.set_index(0, inplace=True)
    df2 = df2.round(5)
    df2.to_csv('data/output/predicted_distances_dataset'+str(dataset)+'.csv')
