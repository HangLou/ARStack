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

    if dataset == 3:
        csv = pd.read_csv(
            'data/input/processed/MSA_Dataset'+str(1)+'.csv')

        remove = ['1X6A', '2XW1', '5M6N', '2A9J', '4ISG', '2GS2',
                  '6UEI', '5PCU', '6PDM', '3T5I', '1S3W', '2QE4']
        keep = set(csv['0'])-set(remove)
        csv = csv.loc[csv['0'].isin(keep)]
        test = pd.read_csv(
            'data/input/processed/MSA_Dataset'+str(dataset)+'.csv')
    else:
        csv = pd.read_csv(
            'data/input/processed/MSA_Dataset'+str(dataset)+'.csv')
        test = pd.read_csv(
            'data/input/processed/MSA_Dataset'+str(dataset)+'.csv')

    #del csv['Unnamed: 0']
    return (csv, test)


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


def predict(X, y, Xt, yt):
    #print(X.shape, y.shape)
    #print(Xt.shape, yt.shape)
    # I create one network to predict whether the magnitude is positive or negative (towards or away from the Baker prediction)
    # Below I transdform the target values for the neural network, I do not have any preferences here.
    y_log = y.copy()
    yt_log = yt.copy()

    #print(min(y_log), min(yt_log))
    #print(max(y_log), max(yt_log))
    #print(min(y), max(y))

    global record
    record = []
    proteins = set(Xt['pdb'])
    print(len(proteins))
    win = 0
    coords = []
    ii = 0

    Xtr = X.sample(frac=1, random_state=9438)
    Xtr = Xtr.iloc[:, :-1]

    ytr = y_log.sample(frac=1, random_state=9438)
    ytr = np.array(ytr)
    ytr = ytr.reshape(-1, 1)

    # model predicts the magnitude
    model = []
    model = Sequential()
    model.add(Dense(200, activation='relu', input_shape=(Xtr.shape[1],)))
    model.add(Dropout(0.2))
    model.add(Dense(200, activation='relu'))
    model.add(Dropout(0.2))
    model.add(Dense(100, activation='softmax'))
    # model.add(Dropout(0.2))
    model.add(Dense(1))

    lr, epcs, bsize = get_parameters()
    sgd = SGD(learning_rate=lr)
    model.compile(loss='mae', optimizer=sgd, metrics=['mse'])

    # 3 epochs and 200 batch size with sgd 0.01 gives equal fotting 37/77
    # 4, 100 works well. Use at least 6 epochs with the settings above
    print(Xtr, ytr)
    model.fit(Xtr, ytr, epochs=epcs, verbose=0, batch_size=bsize)

    for prot in list(proteins):  # [:20]:
        #        print(prot)
     # if prot =='6ILD':
     #   print('going in')
        ii += 1
        print('number: ', ii)
        try:
            Xts = Xt.loc[Xt['pdb'] == prot]
            print(len(Xts))
            yts = np.array(yt.loc[Xts.index])

            yts = yts.reshape(-1, 1)
            Xts = Xts.iloc[:, :-1]

            yp = model.predict(Xts)

            # yp=10**(-1*yp)
            # yp=yp-70
            # yp=transform_back(yp)
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


def train_predict(dataset: int):
    record = []
    csv, test_csv = get_data(dataset)
    df = csv
    df = df[1::3]
    df_test = test_csv
    df_test = df_test[1::3]
    X = df.iloc[:, 6:]
    Xt = df_test.iloc[:, 6:]
    print(X.shape)

    X['pdb'] = df['0']
    Xt['pdb'] = df_test['0']
    y = df['2']
    yt = df_test['2']
    coords, record = predict(X, y, Xt, yt)

    df2 = pd.DataFrame(coords)
    df2.set_index(0, inplace=True)
    df2 = df2.round(5)
    df2.to_csv('data/output/predicted_distances_dataset'+str(dataset)+'.csv')


if __name__ == '__main__':
    for dataset in range(1, 4):

        train_predict(dataset)

    global record
    record = []
    csv, test_csv = get_data(dataset=1)
    df = csv
    df = df[1::3]
    df_test = test_csv
    df_test = df_test[1::3]
    X = df.iloc[:, 6:]
    Xt = df_test.iloc[:, 6:]
    print(X.shape)

    X['pdb'] = df['0']
    Xt['pdb'] = df_test['0']
    y = df['2']
    yt = df_test['2']
    coords, record = predict(X, y, Xt, yt)

    df2 = pd.DataFrame(coords)
    df2.set_index(0, inplace=True)
    df2 = df2.round(5)
    df2.to_csv('data/output/predicted_distances_dataset'+str(dataset)+'.csv')
