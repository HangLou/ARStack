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
from scipy.stats import binom_test
import math
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings('ignore')
print('done importing')


def get_data(dataset=1):
    df = pd.read_csv(
        'data/output/predicted_distances_dataset'+str(dataset)+'.csv')
    df.set_index('0', inplace=True)

#    df=df.iloc[:, int(len(df.columns)/2):]

    return df


def get_vectors(dataset=1):
    #dfv = pd.read_csv('combined_vector_info.csv')
    dfv = pd.read_csv(
        'data/input/processed/vector_Dataset'+str(dataset)+'.csv')
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
    tally_baker = []
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
            # error (alpha coordinates)
            dist5 = np.sqrt((xp[0]-zp[0])**2+(xp[1]-zp[1])**2+(xp[2]-zp[2])**2)
            tally_baker.append(dist5)

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

    return tally_pred, tally_alpha, tally_baker


def plot():
    sns.set()
    dataset_rmsd1 = pd.read_csv('data/output/rmsd_1.csv')
    dataset_rmsd2 = pd.read_csv('data/output/rmsd_2.csv')
    dataset_rmsd3 = pd.read_csv('data/output/rmsd_3.csv')
    dataset_rmsd1['dataset'] = 1
    dataset_rmsd2['dataset'] = 2
    dataset_rmsd3['dataset'] = 3
    MW1 = pd.read_csv('data/output/Dataset_1_MW.csv', names=['MW'])
    MW2 = pd.read_csv('data/output/Dataset_2_MW.csv', names=['MW'])
    MW3 = pd.read_csv('data/output/Dataset_3_MW.csv', names=['MW'])
    MW = pd.concat([MW1, MW2, MW3])
    MW = MW.reset_index()
    MW = MW.rename(columns={'index': 'protein'})
    df = pd.concat([dataset_rmsd1, dataset_rmsd2, dataset_rmsd3])
    df['diff'] = df.iloc[:, 2]-df.iloc[:, 1]

    df = df.rename(columns={'diff': '\u0394RMSD'})
    df = df.rename(columns={'Unnamed: 0': 'protein'})
    df = df.rename(columns={'0': 'stacked', '1': 'alpha', '2': 'baker'})
    df = pd.merge(df, MW, how='left',
                  on='protein').drop_duplicates(keep='first')

    plt.figure(figsize=(8, 6))
    ax = sns.violinplot(data=df, x='dataset',
                        y='\u0394RMSD', linewidth=2.5)
    ax.axhline(0, ls='--', alpha=1, c='k')
    plt.tight_layout()
    plt.savefig("data/output/drmsd_violine_plot.png")

    df['temp'] = df.MW <= 10000
    df['Molecular Weight (W)'] = r'$W \leq 10000$'
    df['Molecular Weight (W)'].loc[df.temp == True] = r'$W \leq 10000$'

    df['temp'] = (df.MW > 10000) & (df.MW <= 20000)
    df['Molecular Weight (W)'].loc[df.temp == True] = r'$10000 < W \leq 20000$'

    df['temp'] = df.MW > 20000
    df['Molecular Weight (W)'].loc[df.temp == True] = r'$20000 < W$'

    plt.figure(figsize=(8, 6))

    ax = sns.violinplot(data=df, x='dataset', y='\u0394RMSD',
                        hue='Molecular Weight (W)', linewidth=2.5)
    #handles, labels = ax.get_legend_handles_labels()

    # ax.legend(handles, [r'$W<=10000$', r'10000<W<=20000',
    #         r'$20000<W$'], loc='upper right')
    ax.axhline(0, ls='--', alpha=1, c='k')
    plt.tight_layout()
    plt.savefig("data/output/drmsd_violine_plot_MW.png")

    #df['Alphafold RMSD > 5'] = df.alpha > 5
    # sns.violinplot(data=df, x='dataset', y='\u0394RMSD',
    #              hue='Alphafold RMSD > 5')
    plt.figure(figsize=(8, 6))
    df['Alphafold'] = df.alpha > 1
    df.Alphafold.loc[df.Alphafold == True] = r'$RMSD > 1$'
    df.Alphafold.loc[df.Alphafold == False] = r'$RMSD \leq 1$'
    ax = sns.violinplot(data=df, x='dataset', y='\u0394RMSD',
                        hue='Alphafold', linewidth=2.5)
    ax.axhline(0, ls='--', alpha=1, c='k')
    plt.tight_layout()
    plt.savefig("data/output/drmsd_violine_plot_cohorts.png")


def evaluation(dataset: int):
    stacked = []
    alpha = []
    baker = []
    df = get_data(dataset)

    win = 0
    loss = 0
    support_vector = get_vectors(dataset)
    # input(support_vector.head())
    proteins = list(df.index)

    for prot in proteins:
        # print(prot)

        try:
            tally1, tally2, tally3 = calc_coord(
                df.loc[prot].dropna(), support_vector.loc[prot].dropna())

            t1 = np.array(tally1)
            t2 = np.array(tally2)
            t3 = np.array(tally3)
            t1 = t1**2
            t2 = t2**2
            t3 = t3**2

            t1 = np.average(t1)  # np.sum(t1)/len(t1)
            t2 = np.average(t2)  # np.sum(t2)/len(t2)
            t3 = np.average(t3)  # np.sum(t2)/len(t2)

            t1 = np.sqrt(t1)
            t2 = np.sqrt(t2)
            t3 = np.sqrt(t3)

            t1 = round(t1, 5)
            t2 = round(t2, 5)
            t3 = round(t3, 5)

            #print('aleph:', t2)
            #print('pred:', t1)
            stacked.append(t1)
            alpha.append(t2)
            baker.append(t3)

            if np.sqrt(t1) < np.sqrt(t2):
                # print('win')
                win += 1
            if np.sqrt(t1) > np.sqrt(t2):
                # print('lose')
                loss += 1

#            input('press enter for next')
        except Exception as e:
            print(traceback.format_exc())
            # input(e)

            nd = 3
    df = pd.DataFrame([stacked, alpha, baker]).T

    df['diff'] = df[0]-df[1]

    df.index = proteins
    df.to_csv('data/output/rmsd_'+str(dataset)+'.csv')
    print('#### Evaluation metric on Dataset ' + str(dataset))
    good = df.loc[df['diff'] < 0]
    bad = df.loc[df['diff'] > 0]
    # df.to_csv('two_per_fam_N_analysis.csv')
    print('mean stacked: ', np.round(
        np.mean(good['diff']), 5), np.std(good['diff']))
    print('median stacked: ', np.round(np.median(good['diff']), 5))
    print('mean alpha: ', np.round(
        np.mean(bad['diff']), 5), np.std(bad['diff']))
    print('median alpha: ', np.round(np.median(bad['diff']), 5))

    print('stacked RMSD mean: ', np.mean(good[0]), np.std(good[0]))
    print('stacked RMSD median : ', np.median(good[0]))
    print('alpha RMSD mean: ', np.mean(bad[1]), np.std(bad[1]))
    print('alpha RMSD median : ', np.median(bad[1]))

    # percentages:
    good_pc = 100*abs(good['diff'])/good[0]
    print('mean stacked pc: ', np.round(np.mean(good_pc), 5), np.std(good_pc))
    print('median stacked pc: ', np.round(np.median(good_pc), 5))
    bad_pc = 100*abs(bad['diff'])/bad[1]
    print('mean alpha pc: ', np.round(np.mean(bad_pc), 5), np.std(bad_pc))
    print('median alpha pc: ', np.round(np.median(bad_pc), 5))

    print('binom p-val: ', binom_test(win, (win+loss), alternative='greater'))

    print('mean stack', np.mean(df[0]), np.std(df[0]))
    print('median stack', np.median(df[0]))
    print('mean alpha', np.mean(df[1]), np.std(df[1]))
    print('median alpha', np.median(df[1]))


if __name__ == '__main__':
    # performance evaluation and save the RMSD of each protein into csv files
    for dataset in range(1, 4):
        evaluation(dataset)
    # visualization
    plot()
