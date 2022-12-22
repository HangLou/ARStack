from collections import Counter
import os#, joblib

def into_df(data):
    global dc
    
    if len(data)>1:
        for n in range(len(data)-1):
            if str(n) in dc:
                cur=dc.get(str(n))
                new=cur+data[n]
                dc[str(n)]=new
            else:
                dc[str(n)]=data[n]


def read_lines(target, folder):
    file = open(target+folder+'/'+folder[:4]+'/msas/uniref90_hits.sto')
    f=file.readlines()
    total=[]
    for line in f:
        if len(line) > 2:
            if "#" not in line:
                total.append(line.strip())

    return (total)


def freq_MSA(total):
    global dc
    dc={}
    count=[]
    new=[]
    s=0
    go=0
    
    for line in total: 
        if 'Chain' in line and '#' not in line:#s!= 0:
            go=1
            positions =[int(pos) for pos, char in enumerate(line.split()[1]) if char != '-']
            into_df(new)
            new=[]

        if '//' in line:
            into_df(new)
            new=[]

        line1 = line.split(' ')[-1]
        line1= list(line1)
        newline=[]
        if go ==1 and '//' not in line:
            for n in positions:
                newline.append(line1[n])
            newline=''.join(newline)
            new.append(newline)
            s=1

    prot_len= len(list(dc.values())[0])
    ls={}
    for i in range(prot_len):
        pos=[]
        for seq in dc.values():
            pos.append(seq[i])
        count=Counter(pos)
        ls[str(i)] = count 

    return(ls)


if __name__=='__main__':
    global df
    target ='../../x_phire/protein_fold_output_round_4/'
#    target ='protein_fold_output_full\\'

    for folder in os.listdir(target):
     try:
        if 'alphafold' in folder:
            total = read_lines(target, folder)
            
            pdb = folder[:4]
            
            ls = freq_MSA(total)

            out =open('MSAs_uniprot_train_set/%s_dict.txt' %pdb, 'w')
            for line in ls.values():
                out.write(str(line))
            out.close()
     except:
         bs=0

#            with open('MSAs/%s_dict.pkl' %pdb, 'wb') as f:
#                joblib.dump(ls, f)

##    dics=[]
##    with open('%s_dict.pkl' %(pdb), 'rb') as f:
##        while True:
##            try:
##                dics.append(joblib.load(f))
##            except EOFError:
##                    break
