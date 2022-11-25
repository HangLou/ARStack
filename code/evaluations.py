from sklearn.metrics import mean_squared_error as mse
import traceback
from Bio.SVDSuperimposer import SVDSuperimposer
import numpy as np
from numpy import array, dot, set_printoptions
from Bio.PDB import *
from Bio import pairwise2
import os
from itertools import groupby
from operator import itemgetter
import warnings
from collections import Counter
import pandas as pd
import math
warnings.filterwarnings('ignore')
print('done importing')


def get_data():
    df = pd.read_csv('train_uniprot_predicted_distances.csv')
    df.set_index('0', inplace=True)

#    df=df.iloc[:, int(len(df.columns)/2):]

    return df


def get_vectors():
    #dfv = pd.read_csv('combined_vector_info.csv')
    dfv = pd.read_csv('JK_two_per_fam_train_vector_info.csv')
    dfv.set_index('0', inplace=True)
    return dfv


def calc_coord(data, vectors):  # called from superimpose
    # data=data*5
    distances = list(data)
    # I do this since I am only looking at one atom from each amino acid at the moment.
    vectors = vectors[1::3]
    # print(len(vectors))
    # print(len(data))
    tally_pred = []
    tally_alpha = []
    for i in range(len(distances)):
        d1 = distances[i]
        values = str(vectors.iloc[i]).strip(']').strip('[').split(',')

        values = [float(i) for i in values]

        direction = np.array(values[6:])  # baker - alphafold values
        zp = np.array(values[:3])  # experimental values
        yp = np.array(values[3:6])  # alohafold values
        dist2 = np.sqrt((yp[0]-zp[0])**2+(yp[1]-zp[1])**2 +
                        (yp[2]-zp[2])**2)  # error (alpha coordinates)
        tally_alpha.append(dist2)

        try:
            xp = direction + yp

    #        angle = math.degrees(math.acos(value)) degree of inverse cos of a value
    #        in a 3D vector, calc using each component (coordinate) separate in relation to the magnitude (length).
            alpha_baker_mag = np.sqrt(
                direction[0]**2+direction[1]**2+direction[2]**2)

            # need angle y first to solve for angle x and y, both methods tried give the same answers...
            angle_y = math.acos(direction[1]/alpha_baker_mag)

            AB = alpha_baker_mag*math.sin(angle_y)

            # these are in radians, not degrees
            angle_x = math.acos(direction[0]/AB)

            angle_z = math.asin(direction[2]/AB)

            # Below calculations are simply using the same angle to find out the coordinates (components) of the predicted distance vector

            new_y = yp[1]+math.cos(angle_y)*d1  # +yp[1]
            AB2 = math.sin(angle_y)*d1

            # +yp[0] #alpha coordinate added to remove us from origo.
            new_x = yp[0]+math.cos(angle_x)*AB2

            new_z = yp[2]+math.sin(angle_z)*AB2  # +yp[2]

            new_coord = [new_x, new_y, new_z]
            # new_coord=xp
            # error (predicted coords)
            dist1 = np.sqrt(
                (new_coord[0]-zp[0])**2+(new_coord[1]-zp[1])**2+(new_coord[2]-zp[2])**2)
            tally_pred.append(dist1)

            # input('halt')

        except Exception as e:
            print(traceback.format_exc())
            # input(e)

            # input(new_coord)

    return tally_pred, tally_alpha


stacked = []
alpha = []


if __name__ == '__main__':

    df = get_data()

    win = 0
    loss = 0
    support_vector = get_vectors()
    # input(support_vector.head())
    proteins = list(df.index)

    for prot in proteins:
        # print(prot)

        try:
            tally1, tally2 = calc_coord(
                df.loc[prot].dropna(), support_vector.loc[prot].dropna())

            t1 = np.array(tally1)
            t2 = np.array(tally2)
            t1 = t1**2
            t2 = t2**2

            t1 = np.average(t1)  # np.sum(t1)/len(t1)
            t2 = np.average(t2)  # np.sum(t2)/len(t2)

            #print('pred:', np.sqrt(t1))
            #print('aleph:', np.sqrt(t2))

            t1 = round(t1, 4)
            t2 = round(t2, 4)

            stacked.append(np.sqrt(t1))
            alpha.append(np.sqrt(t2))

            if np.sqrt(t1) < np.sqrt(t2):
                # print('win')
                win += 1
            if np.sqrt(t1) > np.sqrt(t2):
                # print('lose')
                loss += 1

#            input('press enter for next')
        except Exception as e:
            print(traceback.format_exc())
            input(e)

            nd = 3


print('loss:', loss, '\nwin: ', win)
print('stacked_loss:', np.mean(stacked), 'alpha_loss:', np.mean(alpha))
# input()
