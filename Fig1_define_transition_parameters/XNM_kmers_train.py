# calculate kmer and train markov model for XNM

from hmseekr.kmers import kmers
from hmseekr.train import train

import pandas as pd
import os

# make folder
os.makedirs('./counts', exist_ok=True)
os.makedirs('./markovModels', exist_ok=True)

xheaders = ['XISTint1','XISTrA','XISTrF','XISTint2','XISTrB1','XISTint3','XISTrB2',
            'XISTint4','XISTint5','XISTrD','XISTint6','XISTint7','XISTint8','XISTint9',
            'XISTrE','XISTint10','XISTint11','XISTint12','XISTint13','XISTint14']
nheaders = [f'NEAT1chk{i}' for i in range(1, 23)]
mheaders = [f'MALAT1chk{i}' for i in range(1, 10)]
others = ['XIST-PfIMI', 'MALAT1E', 'MALAT1M', 'ncRNA-a7', 'HOTTIP-MLL', 'FIRRE-RRD',
          'MEG3-NRE','NXF1-enChr','JPX-9','PVT1-22','NKILA-SRSFNRE','e-Ccnd1',
          'e-YY1','e-Med13l','e-Klf6','e-Sp3','e-ID1','PRRX2-eRNA','SEMA3C-eRNA',
          'Blustr-5ss','ASAR6-141chk1','ASAR6-141chk2','ASAR6-141chk3','ASAR6-141chk4','ASAR6-141chk5',
          'HSat3_sense', 'HSat3_antisense']

xnmfiles = xheaders + nheaders + mheaders + others


nulldict = kmers(fadir='./v47.lncRNA.can.500.nodup.OT1fixed.fa',kvec='4',
                 alphabet='ATCG',outputname='v47TSC',
                 outputdir='./counts/')

# read in the optimized transition prob from optimized_transit_p_XNM.csv as a dataframe
transp=pd.read_csv('optimized_transit_p_allfunck4.csv',header=0)


print('finish null dict')

for xnm in xnmfiles:
    print(xnm)
    # count the kmers
    testdict = kmers(fadir=f'./FastaFiles/{xnm}.fa',kvec='4',
                     alphabet='ATCG',outputname=xnm,
                     outputdir='./counts/')
    
    # get the transition prob for the model
    # locate the row index by matching transp col 'file' with xnm
    currentqT = transp.loc[transp['file'] == xnm, 'qT'].iloc[0]
    currentnT = transp.loc[transp['file'] == xnm, 'nT'].iloc[0]

    
    # train the model

    testmodel = train(querydir=f'./counts/{xnm}.dict', nulldir='./counts/v47TSC.dict',
                      kvec='4', alphabet='ATCG', queryT=currentqT, nullT=currentnT,
                      queryPrefix=xnm, nullPrefix='v47TSC', outputdir='./markovModels/')
    
print('finish training')

##################K5


from hmseekr.kmers import kmers
from hmseekr.train import train

import pandas as pd
import os

# make folder
os.makedirs('./k5counts', exist_ok=True)
os.makedirs('./k5markovModels', exist_ok=True)

xheaders = ['XISTint1','XISTrA','XISTrF','XISTint2','XISTrB1','XISTint3','XISTrB2',
            'XISTint4','XISTint5','XISTrD','XISTint6','XISTint7','XISTint8','XISTint9',
            'XISTrE','XISTint10','XISTint11','XISTint12','XISTint13','XISTint14']
nheaders = [f'NEAT1chk{i}' for i in range(1, 23)]
mheaders = [f'MALAT1chk{i}' for i in range(1, 10)]
others = ['XIST-PfIMI', 'MALAT1E', 'MALAT1M', 'ncRNA-a7', 'HOTTIP-MLL', 'FIRRE-RRD',
          'MEG3-NRE','NXF1-enChr','JPX-9','PVT1-22','NKILA-SRSFNRE','e-Ccnd1',
          'e-YY1','e-Med13l','e-Klf6','e-Sp3','e-ID1','PRRX2-eRNA','SEMA3C-eRNA',
          'Blustr-5ss','ASAR6-141chk1','ASAR6-141chk2','ASAR6-141chk3','ASAR6-141chk4','ASAR6-141chk5',
          'HSat3_sense', 'HSat3_antisense']

xnmfiles = xheaders + nheaders + mheaders + others


nulldict = kmers(fadir='./v47.lncRNA.can.500.nodup.OT1fixed.fa',kvec='5',
                 alphabet='ATCG',outputname='v47TSC',
                 outputdir='./k5counts/')

# read in the optimized transition prob from optimized_transit_p_XNM.csv as a dataframe
transp=pd.read_csv('optimized_transit_p_allfunck5.csv',header=0)


print('finish null dict')

for xnm in xnmfiles:
    print(xnm)
    # count the kmers
    testdict = kmers(fadir=f'./fastaFiles/{xnm}.fa',kvec='5',
                     alphabet='ATCG',outputname=xnm,
                     outputdir='./k5counts/')
    
    # get the transition prob for the model
    # locate the row index by matching transp col 'file' with xnm
    currentqT = transp.loc[transp['file'] == xnm, 'qT'].iloc[0]
    currentnT = transp.loc[transp['file'] == xnm, 'nT'].iloc[0]

    
    # train the model

    testmodel = train(querydir=f'./k5counts/{xnm}.dict', nulldir='./k5counts/v47TSC.dict',
                      kvec='5', alphabet='ATCG', queryT=currentqT, nullT=currentnT,
                      queryPrefix=xnm, nullPrefix='v47TSC', outputdir='./k5markovModels/')
    
print('finish training')


##################K6


from hmseekr.kmers import kmers
from hmseekr.train import train

import pandas as pd
import os

# make folder
os.makedirs('./k6counts', exist_ok=True)
os.makedirs('./k6markovModels', exist_ok=True)

xheaders = ['XISTint1','XISTrA','XISTrF','XISTint2','XISTrB1','XISTint3','XISTrB2',
            'XISTint4','XISTint5','XISTrD','XISTint6','XISTint7','XISTint8','XISTint9',
            'XISTrE','XISTint10','XISTint11','XISTint12','XISTint13','XISTint14']
nheaders = [f'NEAT1chk{i}' for i in range(1, 23)]
mheaders = [f'MALAT1chk{i}' for i in range(1, 10)]
others = ['XIST-PfIMI', 'MALAT1E', 'MALAT1M', 'ncRNA-a7', 'HOTTIP-MLL', 'FIRRE-RRD',
          'MEG3-NRE','NXF1-enChr','JPX-9','PVT1-22','NKILA-SRSFNRE','e-Ccnd1',
          'e-YY1','e-Med13l','e-Klf6','e-Sp3','e-ID1','PRRX2-eRNA','SEMA3C-eRNA',
          'Blustr-5ss','ASAR6-141chk1','ASAR6-141chk2','ASAR6-141chk3','ASAR6-141chk4','ASAR6-141chk5',
          'HSat3_sense', 'HSat3_antisense']

xnmfiles = xheaders + nheaders + mheaders + others


nulldict = kmers(fadir='./v47.lncRNA.can.500.nodup.OT1fixed.fa',kvec='6',
                 alphabet='ATCG',outputname='v47TSC',
                 outputdir='./k6counts/')

# read in the optimized transition prob from optimized_transit_p_XNM.csv as a dataframe
transp=pd.read_csv('optimized_transit_p_allfunck6.csv',header=0)


print('finish null dict')

for xnm in xnmfiles:
    print(xnm)
    # count the kmers
    testdict = kmers(fadir=f'./fastaFiles/{xnm}.fa',kvec='6',
                     alphabet='ATCG',outputname=xnm,
                     outputdir='./k6counts/')
    
    # get the transition prob for the model
    # locate the row index by matching transp col 'file' with xnm
    currentqT = transp.loc[transp['file'] == xnm, 'qT'].iloc[0]
    currentnT = transp.loc[transp['file'] == xnm, 'nT'].iloc[0]

    
    # train the model

    testmodel = train(querydir=f'./k6counts/{xnm}.dict', nulldir='./k6counts/v47TSC.dict',
                      kvec='6', alphabet='ATCG', queryT=currentqT, nullT=currentnT,
                      queryPrefix=xnm, nullPrefix='v47TSC', outputdir='./k6markovModels/')
    
print('finish training')

