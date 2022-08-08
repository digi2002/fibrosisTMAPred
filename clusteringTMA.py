from util import scanpycluster_TMAmulti
from plot_util import plotanno,plotsingle,umapplot,plotannotmp
import scanpy as sc
from anno_lungTMA import getOnsampeldata,readMultianno,checknodeanno,genadataAnno
from preprocessTMA import getOnsampleAdata
from util import adatacluster
from lungPara import getPara,getDir,getAnnofile
from util import concatenateadata,getInitalAdata
import scanorama
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def initialCluster(datadir,RES,samplefiles,modeldir,matrixremove = False):
    for fildid in samplefiles.keys():
        modelname = fildid + '_R' + str(RES) + '.h5ad'
        infile = datadir + samplefiles[fildid]
        scanpycluster_TMAmulti(infile, fildid, modeldir + modelname, RES, matrixremove)


def clusterOnsample(modeldir,annofile,fildid,samplefiles,keystrSet,inmodel,RES):
    print('read orignial annotations')
    labelDict = readMultianno(annofile, fildid.replace('lung',''))
    adataonsample = getOnsampleAdata(datadir, fildid, samplefiles, inmodel)
    adata = genadataAnno(labelDict, adataonsample, fildid)
    modelname = modeldir + fildid + '_R' + str(RES) + '_onsample.h5ad'
    ad = sc.read(modelname)
    for keystr in keystrSet.keys():
        plotanno(ad, keystr, figdir, fildid, RES,'_onsample'+keystr)
    plotanno(ad, 'label', figdir, fildid, RES, '_onsample' + 'label')


def clustersample(modeldir,annofile,fildid,samplefiles,keystrSet,inmodel,RES):
    print('read orignial annotations')
    labelDict = readMultianno(annofile, fildid.replace('lung',''))
    adatasample = sc.read(inmodel)
    adata = genadataAnno(labelDict, adatasample, fildid)
    modelname = modeldir + fildid + '_R' + str(RES) + '_allpixel.h5ad'
    adatacluster(adata, RES, modelname)
    ad = sc.read(modelname)
    for keystr in keystrSet.keys():
        plotanno(ad, keystr, figdir, fildid, RES,'_allpixel'+keystr)
    plotanno(ad, 'label', figdir, fildid, RES, '_allpixel' + 'label')


def clusterallOnsample(modeldir,annofile,samplefiles,RES,batch=False):
    adatas=[]
    for fildid in samplefiles.keys():
        print('read orignial annotations')
        labelDict = readMultianno(annofile, fildid.replace('lung', ''))
        inmodel = modeldir + fildid + '_R' + str(0.1) + '.h5ad'
        adataonsample = getOnsampleAdata(datadir, fildid, samplefiles, inmodel)
        adata = genadataAnno(labelDict, adataonsample, fildid)
        adatas.append(adata)
    if batch==True:
        adatas = scanorama.correct_scanpy(adatas, return_dimred=True)
        alladata = concatenateadata(adatas)
        modelname = modeldir + 'alllung' + '_R' + str(RES) + '_onsample_batchcorrect.h5ad'
    else:
        alladata = concatenateadata(adatas)
        modelname = modeldir + 'alllung' + '_R' + str(RES) + '_onsample.h5ad'
    adatacluster(alladata, RES, modelname,batch)

def getTMAOnsample(modeldir,datadir,samplefiles,annofile,fildid,RES):
    print('read orignial annotations')
    labelDict = readMultianno(annofile, fildid.replace('lung', ''))
    inmodel = modeldir + fildid + '_R' + str(RES) + '.h5ad'
    adataonsample = getOnsampleAdata(datadir, fildid, samplefiles, inmodel)

    return adataonsample

def getAdataanno(datadir,samplefiles,annofile,fildid,RES):
    print('read orignial annotations')
    labelDict = readMultianno(annofile, fildid.replace('lung', ''))
    infile = datadir + samplefiles[fildid]
    adata = getInitalAdata(infile, fildid, num=0)
    adata = genadataAnno(labelDict, adata, fildid)
    return adata

def scatterUmapcore(df):
    coord=np.array(df['umap'].tolist())
    color=df['leiden'].tolist()
    color=[int(i) for i in color]
    plt.scatter(coord[:,0],coord[:,1],c=color)
    plt.show()


def plotCoreumap(ad):
    umapcoordinates = ad.obsm['X_umap'].tolist()
    leidens = ad.obs['leiden'].tolist()
    labels = ad.obs['label'].tolist()
    df = pd.DataFrame({'umap': umapcoordinates, 'leiden': leidens, 'labels': labels})
    for sellabel in set(labels):
        newdf = df.loc[df['labels'] == sellabel]
        scatterUmapcore(newdf)

def plotallonsample(modeldir,RES,batch=False):
    if batch==True:
        modelname = modeldir + 'alllung' + '_R' + str(RES) + '_onsample_batchcorrect.h5ad'
    else:
        modelname = modeldir + 'alllung' + '_R' + str(RES) + '_onsample.h5ad'
    ad = sc.read(modelname)
    umapplot(ad, 'leiden')
    plotCoreumap(ad)


def labeladata(annofile,fildid,keystrSet,inmodel,RES):
    print('read orignial annotations')
    labelDict = readMultianno(annofile, fildid.replace('lung',''))
    adatasample = sc.read(inmodel)
    adata = genadataAnno(labelDict, adatasample, fildid)
    adata.write(inmodel)
    for keystr in keystrSet.keys():
        plotannotmp(adata, keystr, figdir, fildid, RES, '_allpixel'+keystr)
    plotannotmp(adata, 'label', figdir, fildid, RES, '_allpixel' + 'label')


def scatterOnsamplecore(ad):
    spatialcoordinates = ad.obsm['spatial'].tolist()
    labels = ad.obs['label'].tolist()
    ages = ad.obs['age'].tolist()
    df = pd.DataFrame({'spatial': spatialcoordinates, 'labels': labels,'age':ages})
    coord=np.array(df['spatial'].tolist())
    color=df['age'].tolist()
    color=[int(i) for i in color]
    plt.scatter(coord[:,0],(-1)*coord[:,1],c=color,s = 2)
    plt.show()

if __name__ == "__main__":
    datadir, modeldir, figdir, resdir, metadatadir, coredatadir = getDir()
    samplefiles, keystrSet, RES = getPara()
    annofile = getAnnofile()
    #keystr='histology'
    removematrix =False
    app = ''

    initialCluster(datadir,RES,samplefiles, modeldir,matrixremove = removematrix)

    for sampleid in samplefiles.keys():
        modelname = sampleid + '_R' + str(RES) + '.h5ad'
        adata = sc.read(modeldir +modelname)
        plotsingle(adata,sampleid,RES,figdir,'leiden',app= app)
        plotsingle(adata, sampleid, RES, figdir, 'kmeans2', app=app)
        plotsingle(adata, sampleid, RES, figdir, 'kmeans3', app=app)
        plotsingle(adata, sampleid, RES, figdir, 'kmeans4', app=app)
        plotsingle(adata, sampleid, RES, figdir, 'kmeans5', app=app)
        plotsingle(adata, sampleid, RES, figdir, 'kmeans6', app=app)
        plotsingle(adata, sampleid, RES, figdir, 'kmeans7', app=app)
        plotsingle(adata, sampleid, RES, figdir, 'kmeans8', app=app)
        plotsingle(adata, sampleid, RES, figdir, 'kmeans9', app=app)
        plotsingle(adata, sampleid, RES, figdir, 'kmeans10', app=app)





