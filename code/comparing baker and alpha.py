import numpy as np
import pandas as pd
import traceback


def get_vectors():
    #dfv=pd.read_csv('JK_test_vector_info.csv')
    dfv=pd.read_csv('JK_test_vector_info_upicif.csv')
    dfvp=pd.read_csv('JK_test_2_vector_info_upicif.csv')
    
    dfv.set_index('0', inplace=True)
    dfvp.set_index('0', inplace=True)

    dfv=pd.concat([dfv,dfvp],axis=0,sort=False)

    good=['7QUE', '7TBS', '7TM8', '7T8O', '7QPU', '7WH0', '7WGJ', '7TRW', '5SD6', '7T2Z', '7TWE', '7TMV', '7QQ6', '7WKQ', '7WZY', '7QU0', '7TAV', '7TM4', '7TCM', '7TMU', '7QU3', '7QUW', '7TB5', '7WH4', '7TFM', '7TM9', '7TT9', '7TXS', '7W74', '7TB7', '5SDQ', '7TLR', '7WBR', '7WRP', '7TVD', '7TYE', '7QY3', '7T79', '7T7S', '7WH1', '7TQ8', '7U5Q', '7U1H', '7TMD', '7YXE', '7YWR', '7W83', '7R4N', '7QNO', '7TA9', '7WKC', '7T4X', '7U56', '7TG5', '7T36', '7T3Y', '7TMG', '7TZP', '7WSV', '7U5F', '7T35', '7TA5', '7TE0', '7QOS', '7QPZ', '7WA8', '5SCS', '7T29', '7U35', '7QYI', '7TI7', '7QJK', '7QRO', '7W6G', '7TE7', '7U5Y', '7WAB', '7T88', '7WF5', '7THI', '7TO6', '7U28', '7T1Q', '7TM7', '7T7J', '7T5Y', '7WKH', '7QWM', '7TGQ', '7WSJ', '7WWH', '7WAN', '7W6F', '7T71', '7TWZ', '7TKV', '7QGA', '7TOC', '7TBG', '7TMB', '7T39', '7YZV', '7TEM', '7WBD', '7QFE', '7QHG', '7QGJ', '7W7J', '7QGP', '7WIN', '7R5H', '7T4L']

    dfv=dfv.loc[good]

    return dfv

def calc_coord(vectors): #called from superimpose
    vectors = list(vectors)[1:]#[1::3] # I do this since I am only looking at one atom from each amino acid at the moment.
    tally_alpha=[]
    tally_baker=[]
    for i in range(len(vectors)):
        values = str(vectors[i]).strip(']').strip('[').split(',')
        values = [float(i) for i in values]
        direction = np.array(values[6:]) # baker - alphafold values
        zp=np.array(values[:3]) # experimental values

        yp=np.array(values[3:6]) # alphafold values
        dist2 = np.sqrt((yp[0]-zp[0])**2+(yp[1]-zp[1])**2+(yp[2]-zp[2])**2) #error (alpha coordinates)
        tally_alpha.append(dist2)

        xp = direction + yp
        dist5 = np.sqrt((xp[0]-zp[0])**2+(xp[1]-zp[1])**2+(xp[2]-zp[2])**2) #error (alpha coordinates)
        tally_baker.append(dist5)

    return tally_alpha, tally_baker


alpha=[]
baker=[]

if __name__ == '__main__':

    win=0
    loss=0

    support_vector = get_vectors()

    proteins = list(support_vector.index)

    for prot in proteins:
        try:
            tally1, tally2 = calc_coord(support_vector.loc[prot].dropna())
            t1=np.array(tally1)
            t2=np.array(tally2)

            t1=t1**2
            t2=t2**2

            t1=np.average(t1)#np.sum(t1)/len(t1)
            t2=np.average(t2)#np.sum(t2)/len(t2)

            t1=np.sqrt(t1)
            t2=np.sqrt(t2)

            t1=round(t1,5)
            t2=round(t2,5)

            alpha.append(t1)
            baker.append(t2)

            if t1 > t2:
                win+=1

            else:
                loss+=1

        except Exception as e:
            print(traceback.format_exc())
            input(e)

print('loss:', loss, '\nwin: ', win)
