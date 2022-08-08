import scanpy as sc
import numpy as np
from collections import Counter


dir = './'
modeldir='../maldiclustering/models_mice/'
print(modeldir)
figdir=dir+'figures_fibrosisTMA/'
batch=False

fibrosis = {'COVID_A42': ['2','5'], 'FibrosisC_A4': ['13','0','2'],'FibrosisB_A6':['0']}
nonfibrosis={'COVID_A42': ['1','4','6','8','9','10','13'],
             'FibrosisC_A4': ['3','6','7','19'],
             'FibrosisB_A6':['2','4','11'],
             'Flu_A20':['0','3'],'Flu_A19':['10','9']}
features = ['1809.6495', '1444.5138', '1791.6316', '2122.7349', '2057.7407', '1793.6561']

testset={'COVID_A43':['0','1','3','4','6','8','12'],
         'FibrosisB_B5':['3','10','0'],
         'FibrosisB_A9':['1'],
         'FibrosisB_A13':['2']
         }

testsetname={'COVID_A43':
                 {'0':'1_fibrosis',
                  '1':'1_fibrosis',
                  '3':'1_fibrosis',
                  '4':'1_fibrosis',
                  '6':'2_nonfib',
                  '8':'2_nonfib',
                  '12':'2_nonfib'},
             'FibrosisB_B5':{
                 '3':'2_nonfib',
                 '10':'2_nonfib',
                 '0':'1_fibrosis'},
            'FibrosisB_A9':{'1':'1_fibrosis'},
            'FibrosisB_A13':{'2':'1_fibrosis'}
            }




def genTrainFibronew(fildid,features):
    #fildid 是normal tissue
    #fibrosis 1
    #non-fibrosis 2
    #normal 0
    np.random.seed(0)
    X=[]
    y=[]
    X1=[]
    y1=[]
    featlist=[]

    #制作fibrosis sample
    for fib in fibrosis.keys():
        adata = sc.read(modeldir + fib + 'DE.h5ad')
        featlist=adata.var_names.tolist()
        ad = adata[adata.obs['leiden'].isin(fibrosis[fib]), :].copy()
        if len(features)!=0:
            ad = ad[:,ad.var_names.isin(features)].copy()
        X.extend(np.array(ad.X))
    y=list(np.ones(len(X)))
    positiveNum=len(X)

    if len(nonfibrosis)!=0:
        for fib in nonfibrosis.keys():
            adata = sc.read(modeldir + fib + 'DE.h5ad')
            ad = adata[adata.obs['leiden'].isin(nonfibrosis[fib]), :].copy()
            if len(features) != 0:
                ad = ad[:,ad.var_names.isin(features)].copy()
            X1.extend(np.array(ad.X))
        X1=np.array(X1)
        if fildid!='':
            #indices = np.random.choice(X1.shape[0], int(positiveNum/2), replace=False)
            indices = np.random.choice(X1.shape[0], int(positiveNum/2), replace=False)
        elif fildid=='':
            indices = np.random.choice(X1.shape[0], int(positiveNum), replace=False)
        X1 = X1[indices]
        for i in range(len(X1)):
            y1.append(2)
        X.extend(X1)
        y.extend(y1)


    if fildid!='':
        adata = sc.read(modeldir + fildid + 'DE.h5ad')
        if len(features) != 0:
            adata = adata[:, adata.var_names.isin(features)].copy()
        print(adata)
        normal=np.array(adata.X)
        if len(nonfibrosis)!=0:
            indices = np.random.choice(normal.shape[0], int(positiveNum/2), replace=False)
            #indices = np.random.choice(normal.shape[0], int(positiveNum/2), replace=False)
        else:
            indices = np.random.choice(normal.shape[0], int(positiveNum), replace=False)
        normal=normal[indices]
        X.extend(normal)
        y.extend(np.zeros(len(normal)))
    print(Counter(y))
    return X,y,featlist