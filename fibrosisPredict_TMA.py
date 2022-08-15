#©2022 Qi Sun


import numpy as np
from sklearn.svm import NuSVC,LinearSVC
import pickle,os
from sklearn.preprocessing import MinMaxScaler,StandardScaler

from trainDataset import genTrainFibronew
from collections import Counter
from tmaPara import getPara,getDir,getAnnofile
from DEfeature import genTMAfeat
from plot_util import plotsingle,singlespatial
from clusteringTMA import getTMAOnsample,getAdataanno
from anno_lungTMA import readMultianno,genadataAnno
from preprocessTMA import getOnsampleAdata

def classifyROC(adata,modelname,features,rmmatrix=False):
    scaler = StandardScaler()
    normalfildid = ''
    X_train, y_train,featlist = genTrainFibronew(normalfildid, features)
    print(Counter(y_train))
    scaler.fit(X_train)
    X_train = scaler.transform(X_train)
    print('training '+ modelname)
    clf = NuSVC(random_state=0,nu=0.01)
    clf.fit(X_train, y_train)
    pickle.dump(clf, open(modeldir+modelname, 'wb'))
    clf = pickle.load(open(modeldir+modelname, 'rb'))
    if len(features)!=0:
        adata = adata[:, adata.var_names.isin(features)].copy()
    X_test = list(np.array(adata.X))
    X_test = scaler.transform(X_test)
    test_preds = clf.predict(X_test)
    print(Counter(test_preds))
    labels = []
    for i in test_preds:
        if i == 0:
            labels.append('0_normal')
        elif i == 1:
            labels.append('1_fibrosis')
        elif i == 2:
            labels.append('2_nonfib')
    adata.obs['class'] = labels
    print(Counter(labels))
    if rmmatrix == True:
        app = 'nonormal_rmmatrix_top_'+str(len(features))
    else:
        app = 'nonormal_top_' + str(len(features))

    singlespatial(adata, '1B', RES, figdir, 'class', str(Counter(labels)),app = app)


if __name__ == "__main__":
    datadir, modeldir, figdir, resdir, metadatadir, coredatadir = getDir()
    samplefiles, keystrSet, RES = getPara()
    annofile = getAnnofile()
    fildid = 'lung1B'
    modelname = 'test'
    rmmatrix = False

    for topN in [3]:
        if rmmatrix:
            labelDict = readMultianno(annofile, fildid.replace('lung', ''))
            inmodel = modeldir + fildid + '_R' + str(RES) + '.h5ad'
            adataonsample = getOnsampleAdata(datadir, fildid, samplefiles, inmodel)
            ad = genadataAnno(labelDict, adataonsample, fildid)
        else:
            ad = getAdataanno(datadir,samplefiles,annofile,fildid,RES)
        tops, adata = genTMAfeat(topN,ad)
        classifyROC(ad, modelname,tops,rmmatrix)
        print(tops)

