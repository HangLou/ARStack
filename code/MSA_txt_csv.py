import numpy as np
import pandas as pd
import os

numbers = np.arange(0, 1000, 100)
numbers = dict(zip(numbers, numbers))


#df_unique=pd.read_csv('test_unique_proteins.csv', header = None)
# ls=df_unique[1].values


ko = 0
for file in os.listdir('data/input/raw/MSAs/'):
 # if file[:4] in ls:
    ko += 1
    if ko in numbers:
        print(ko)
    maj_dc = {}
    f = open('data/input/raw/MSAs/'+file, 'r').readlines()
    pdb = file[:4]
    for line in f:
        line = line.split('Counter')
        line = line[1:]
        for n in range(len(line)):
            dc = {}
            newline = line[n].split('{')[1].split('}')[0].split(',')
            for m in newline:
                m = m.split(':')
#                input(m[0].strip().strip('"'))
                dc[m[0].strip().strip("'")] = int(m[1])
                maj_dc[n] = dc

    df = pd.DataFrame(maj_dc).T
    df.fillna(0, inplace=True)
    thresh = np.sum(df.iloc[:, 0])/len(df)
    if 1 == 1:  # thresh > 100:
        # try:
        #    del df['-']
        # except:
        #    bd=0
        df = df.T
        df = df/np.sum(df)
        df = df.T
#        print(np.sum(df.iloc[0,:]))
        # input(df.head())

        df.to_csv('data/input/raw/MSA_csv/%s_MSA_tuples.csv' % pdb)
