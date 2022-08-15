#©2022 Qi Sun


import numpy as np
import scanpy as sc
import pandas as pd
from sklearn.cluster import KMeans
from DEfeature import readmatrix
import scanorama

def concatenateadata(adata_spatial):
    if len(adata_spatial) == 1:
        adata_spatial=adata_spatial[0]
    elif len(adata_spatial)==2:
        adata_spatial = adata_spatial[0].concatenate(
            adata_spatial[1],
            batch_key="library_id",
            uns_merge="unique",
            batch_categories=[
                k
                for d in [
                    adata_spatial[0].uns["spatial"],
                    adata_spatial[1].uns["spatial"]
                ]
                for k, v in d.items()
            ]
        )
    elif len(adata_spatial)==5:
        adata_spatial = adata_spatial[0].concatenate(
            adata_spatial[1],adata_spatial[2], adata_spatial[3], adata_spatial[4],
            batch_key="library_id",
            uns_merge="unique",
            batch_categories=[
                k
                for d in [
                    adata_spatial[0].uns["spatial"],
                    adata_spatial[1].uns["spatial"],
                    adata_spatial[2].uns["spatial"],
                    adata_spatial[3].uns["spatial"],
                    adata_spatial[4].uns["spatial"]
                ]
                for k, v in d.items()
            ]
        )
    elif len(adata_spatial) == 4:
        adata_spatial = adata_spatial[0].concatenate(
            adata_spatial[1], adata_spatial[2], adata_spatial[3],
            batch_key="library_id",
            uns_merge="unique",
            batch_categories=[
                k
                for d in [
                    adata_spatial[0].uns["spatial"],
                    adata_spatial[1].uns["spatial"],
                    adata_spatial[2].uns["spatial"],
                    adata_spatial[3].uns["spatial"]
                ]
                for k, v in d.items()
            ]
        )
    elif len(adata_spatial) == 3:
        adata_spatial = adata_spatial[0].concatenate(
            adata_spatial[1], adata_spatial[2],
            batch_key="library_id",
            uns_merge="unique",
            batch_categories=[
                k
                for d in [
                    adata_spatial[0].uns["spatial"],
                    adata_spatial[1].uns["spatial"],
                    adata_spatial[2].uns["spatial"]
                ]
                for k, v in d.items()
            ]
        )
    elif len(adata_spatial) == 10:
        adata_spatial = adata_spatial[0].concatenate(
            adata_spatial[1], adata_spatial[2],adata_spatial[3], adata_spatial[4],
            adata_spatial[5], adata_spatial[6], adata_spatial[7], adata_spatial[8],
            adata_spatial[9],
            batch_key="library_id",
            uns_merge="unique",
            batch_categories=[
                k
                for d in [
                    adata_spatial[0].uns["spatial"],
                    adata_spatial[1].uns["spatial"],
                    adata_spatial[2].uns["spatial"],
                    adata_spatial[3].uns["spatial"],
                    adata_spatial[4].uns["spatial"],
                    adata_spatial[5].uns["spatial"],
                    adata_spatial[6].uns["spatial"],
                    adata_spatial[7].uns["spatial"],
                    adata_spatial[8].uns["spatial"],
                    adata_spatial[9].uns["spatial"]
                ]
                for k, v in d.items()
            ]
        )
    elif len(adata_spatial) == 9:
        adata_spatial = adata_spatial[0].concatenate(
            adata_spatial[1], adata_spatial[2], adata_spatial[3], adata_spatial[4],
            adata_spatial[5], adata_spatial[6], adata_spatial[7], adata_spatial[8],
            batch_key="library_id",
            uns_merge="unique",
            batch_categories=[
                k
                for d in [
                    adata_spatial[0].uns["spatial"],
                    adata_spatial[1].uns["spatial"],
                    adata_spatial[2].uns["spatial"],
                    adata_spatial[3].uns["spatial"],
                    adata_spatial[4].uns["spatial"],
                    adata_spatial[5].uns["spatial"],
                    adata_spatial[6].uns["spatial"],
                    adata_spatial[7].uns["spatial"],
                    adata_spatial[8].uns["spatial"]
                ]
                for k, v in d.items()
            ]
        )
    elif len(adata_spatial) == 8:
        adata_spatial = adata_spatial[0].concatenate(
            adata_spatial[1], adata_spatial[2], adata_spatial[3], adata_spatial[4],
            adata_spatial[5], adata_spatial[6], adata_spatial[7],
            batch_key="library_id",
            uns_merge="unique",
            batch_categories=[
                k
                for d in [
                    adata_spatial[0].uns["spatial"],
                    adata_spatial[1].uns["spatial"],
                    adata_spatial[2].uns["spatial"],
                    adata_spatial[3].uns["spatial"],
                    adata_spatial[4].uns["spatial"],
                    adata_spatial[5].uns["spatial"],
                    adata_spatial[6].uns["spatial"],
                    adata_spatial[7].uns["spatial"]
                ]
                for k, v in d.items()
            ]
        )
    elif len(adata_spatial)==6:
        adata_spatial = adata_spatial[0].concatenate(
            adata_spatial[1],adata_spatial[2], adata_spatial[3], adata_spatial[4],adata_spatial[5],
            batch_key="library_id",
            uns_merge="unique",
            batch_categories=[
                k
                for d in [
                    adata_spatial[0].uns["spatial"],
                    adata_spatial[1].uns["spatial"],
                    adata_spatial[2].uns["spatial"],
                    adata_spatial[3].uns["spatial"],
                    adata_spatial[4].uns["spatial"],
                    adata_spatial[5].uns["spatial"]
                ]
                for k, v in d.items()
            ]
        )

    elif len(adata_spatial) == 7:
        adata_spatial = adata_spatial[0].concatenate(
            adata_spatial[1], adata_spatial[2], adata_spatial[3], adata_spatial[4],
            adata_spatial[5], adata_spatial[6],
            batch_key="library_id",
            uns_merge="unique",
            batch_categories=[
                k
                for d in [
                    adata_spatial[0].uns["spatial"],
                    adata_spatial[1].uns["spatial"],
                    adata_spatial[2].uns["spatial"],
                    adata_spatial[3].uns["spatial"],
                    adata_spatial[4].uns["spatial"],
                    adata_spatial[5].uns["spatial"],
                    adata_spatial[6].uns["spatial"]
                ]
                for k, v in d.items()
            ]
        )


    elif len(adata_spatial)==15:
        adata_spatial = adata_spatial[0].concatenate(
            adata_spatial[1],
            adata_spatial[2], adata_spatial[3], adata_spatial[4],
            adata_spatial[5], adata_spatial[6], adata_spatial[7], adata_spatial[8],
            adata_spatial[9], adata_spatial[10], adata_spatial[11], adata_spatial[12],
            adata_spatial[13], adata_spatial[14],

            batch_key="library_id",
            uns_merge="unique",
            batch_categories=[
                k
                for d in [
                    adata_spatial[0].uns["spatial"],
                    adata_spatial[1].uns["spatial"],
                    adata_spatial[2].uns["spatial"],
                    adata_spatial[3].uns["spatial"],
                    adata_spatial[4].uns["spatial"],
                    adata_spatial[5].uns["spatial"],
                    adata_spatial[6].uns["spatial"],
                    adata_spatial[7].uns["spatial"],
                    adata_spatial[8].uns["spatial"],
                    adata_spatial[9].uns["spatial"],
                    adata_spatial[10].uns["spatial"],
                    adata_spatial[11].uns["spatial"],
                    adata_spatial[12].uns["spatial"],
                    adata_spatial[13].uns["spatial"],
                    adata_spatial[14].uns["spatial"]
                ]
                for k, v in d.items()
            ],
        )

    elif len(adata_spatial)==22:
        adata_spatial = adata_spatial[0].concatenate(
            adata_spatial[1],
            adata_spatial[2], adata_spatial[3], adata_spatial[4],
            adata_spatial[5], adata_spatial[6], adata_spatial[7], adata_spatial[8],
            adata_spatial[9], adata_spatial[10], adata_spatial[11], adata_spatial[12],
            adata_spatial[13], adata_spatial[14], adata_spatial[15], adata_spatial[16],
            adata_spatial[17], adata_spatial[18], adata_spatial[19], adata_spatial[20],
            adata_spatial[21],

            batch_key="library_id",
            uns_merge="unique",
            batch_categories=[
                k
                for d in [
                    adata_spatial[0].uns["spatial"],
                    adata_spatial[1].uns["spatial"],
                    adata_spatial[2].uns["spatial"],
                    adata_spatial[3].uns["spatial"],
                    adata_spatial[4].uns["spatial"],
                    adata_spatial[5].uns["spatial"],
                    adata_spatial[6].uns["spatial"],
                    adata_spatial[7].uns["spatial"],
                    adata_spatial[8].uns["spatial"],
                    adata_spatial[9].uns["spatial"],
                    adata_spatial[10].uns["spatial"],
                    adata_spatial[11].uns["spatial"],
                    adata_spatial[12].uns["spatial"],
                    adata_spatial[13].uns["spatial"],
                    adata_spatial[14].uns["spatial"],
                    adata_spatial[15].uns["spatial"],
                    adata_spatial[16].uns["spatial"],
                    adata_spatial[17].uns["spatial"],
                    adata_spatial[18].uns["spatial"],
                    adata_spatial[19].uns["spatial"],
                    adata_spatial[20].uns["spatial"],
                    adata_spatial[21].uns["spatial"]
                ]
                for k, v in d.items()
            ],
        )
    return adata_spatial


def removeMatrixfeat(adata):
    matrixlist = readmatrix()
    matrixlist = [float(i) for i in matrixlist]
    matrixlist = [str(i) for i in matrixlist]
    featlist=adata.var_names.tolist()
    newfeatlist = [i for i in featlist if i not in matrixlist]
    adata = adata[:, adata.var_names.isin(newfeatlist)].copy()
    return adata


def getInitalAdata(infile,fildid,num=0):
    xs, ys, nodefeatlabels, nodefeatlabelsdict, nodefeats = readOrifile(infile, num)
    nodefeatlabels = [float(i) for i in nodefeatlabels]
    nodefeatlabels = [str(i) for i in nodefeatlabels]
    varindex = pd.DataFrame(index=nodefeatlabels)
    adata = genadata(xs, ys, nodefeats, fildid, varindex)
    sc.pp.normalize_total(adata, target_sum=1, inplace=True)
    sc.pp.log1p(adata, base=2)
    return adata

def scanpycluster_TMAmulti(infile, fildid,modelname, RES,matrixremove=False):
    num = 0
    xs, ys, nodefeatlabels, nodefeatlabelsdict, nodefeats = readOrifile(infile, num)
    nodefeatlabels = [float(i) for i in nodefeatlabels]
    nodefeatlabels = [str(i) for i in nodefeatlabels]
    varindex = pd.DataFrame(index=nodefeatlabels)
    adata = genadata(xs, ys, nodefeats, fildid, varindex)
    sc.pp.normalize_total(adata, target_sum=1, inplace=True)
    sc.pp.log1p(adata, base=2)

    if matrixremove == True:
        adata = removeMatrixfeat(adata)
    adata.raw = adata
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=RES)
    '''
    kmeans = KMeans(n_clusters=2, random_state=0).fit(adata.X)
    adata.obs['kmeans2'] = kmeans.labels_.astype(str)
    kmeans = KMeans(n_clusters=3, random_state=0).fit(adata.X)
    adata.obs['kmeans3'] = kmeans.labels_.astype(str)
    kmeans = KMeans(n_clusters=4, random_state=0).fit(adata.X)
    adata.obs['kmeans4'] = kmeans.labels_.astype(str)
    kmeans = KMeans(n_clusters=5, random_state=0).fit(adata.X)
    adata.obs['kmeans5'] = kmeans.labels_.astype(str)
    kmeans = KMeans(n_clusters=6, random_state=0).fit(adata.X)
    adata.obs['kmeans6'] = kmeans.labels_.astype(str)
    kmeans = KMeans(n_clusters=7, random_state=0).fit(adata.X)
    adata.obs['kmeans7'] = kmeans.labels_.astype(str)
    kmeans = KMeans(n_clusters=8, random_state=0).fit(adata.X)
    adata.obs['kmeans8'] = kmeans.labels_.astype(str)
    kmeans = KMeans(n_clusters=9, random_state=0).fit(adata.X)
    adata.obs['kmeans9'] = kmeans.labels_.astype(str)
    kmeans = KMeans(n_clusters=10, random_state=0).fit(adata.X)
    adata.obs['kmeans10'] = kmeans.labels_.astype(str)
    '''
    print(adata)
    adata.write(modelname)

def scanpycluster(dir,fildid,samplefiles,modelname,RES):
    num=0
    infile = dir + samplefiles[fildid]
    xs, ys, nodefeatlabels, nodefeatlabelsdict, nodefeats = readOrifile(infile, num)

    varindex = pd.DataFrame(index=nodefeatlabels)
    adata = genadata(xs, ys, nodefeats, fildid, varindex)
    sc.pp.normalize_total(adata, target_sum=1, inplace=True)
    sc.pp.log1p(adata, base=2)
    adata.raw = adata
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata,resolution=RES)
    print(adata)
    adata.write(modelname)

def adatacluster(adata,RES,modelname,batch=False):
    if batch==True:
        sc.pp.neighbors(adata, use_rep="X_scanorama")
    else:
        sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=RES)
    kmeans = KMeans(n_clusters=10, random_state=0).fit(adata.X)
    adata.obs['kmeans10'] = kmeans.labels_.astype(str)
    kmeans = KMeans(n_clusters=5, random_state=0).fit(adata.X)
    adata.obs['kmeans5'] = kmeans.labels_.astype(str)
    print(adata)
    adata.write(modelname)


def genadata(xs,ys,nodefeats,lidid,varindex):
    pos=np.stack((xs,ys), axis=-1)
    adata= sc.AnnData(X=np.array(nodefeats))
    adata.obsm['spatial'] = pos
    adata.uns['spatial'] = {lidid: {}}
    adata.obs['spatial'] = lidid
    adata.var = varindex
    return adata

def readOrifile(filename,num=0):
    nodefeatlabelsdict={}
    cnt=0
    f=open(filename)
    pre=f.readline()
    pre=f.readline()
    pre=f.readline()
    line = f.readline()
    tokens = line.strip().split('\t')
    nodefeatlabels = [str(i) for i in tokens]
    for i in range(len(nodefeatlabels)):
        nodefeatlabelsdict[float(nodefeatlabels[i])]=i
    xs=[]
    ys=[]
    nodefeats=[]
    cnt=0
    for line in f:
        if num!=0:
            cnt=cnt+1
            if cnt>num:
                break
        if line.strip()!='':
            tokens=line.strip().split('\t')
            x=float(tokens[1])
            y=float(tokens[2])
            feat=tokens[3:-2]
            nodefeat = [float(k) for k in feat]
            xs.append(x)
            ys.append(y)
            nodefeats.append(nodefeat)
    return np.array(xs),np.array(ys),nodefeatlabels,nodefeatlabelsdict,np.array(nodefeats)


def scanpyIntecluster_multi(dir,samples,samplefiles,batch,modelname,RES):
    num = 0
    adatas = []
    for sample in samples:
        infile = dir + samplefiles[sample]
        print(infile)
        xs, ys, nodefeatlabels, nodefeatlabelsdict, nodefeats = readOrifile(infile, num)
        varindex = pd.DataFrame(index=nodefeatlabels)
        adata = genadata(xs, ys, nodefeats, sample, varindex)
        sc.pp.normalize_total(adata, target_sum=1, inplace=True)
        sc.pp.log1p(adata, base=2)
        adatas.append(adata)
    if batch == True:
        adata_spatial = scanorama.correct_scanpy(adatas, return_dimred=True)
    else:
        adata_spatial = adatas
    adata_spatial = concatenateadata(adata_spatial)
    adata_spatial.raw = adata_spatial
    if batch == True:
        sc.pp.neighbors(adata_spatial, use_rep="X_scanorama")
    else:
        sc.pp.neighbors(adata_spatial)
    sc.tl.umap(adata_spatial,min_dist=1)
    sc.tl.leiden(adata_spatial, resolution=RES)

    kmeans = KMeans(n_clusters=2, random_state=0).fit(adata_spatial.X)
    adata_spatial.obs['kmeans2'] = kmeans.labels_.astype(str)

    kmeans = KMeans(n_clusters=3, random_state=0).fit(adata_spatial.X)
    adata_spatial.obs['kmeans3'] = kmeans.labels_.astype(str)

    kmeans = KMeans(n_clusters=4, random_state=0).fit(adata_spatial.X)
    adata_spatial.obs['kmeans4'] = kmeans.labels_.astype(str)

    kmeans = KMeans(n_clusters=5, random_state=0).fit(adata_spatial.X)
    adata_spatial.obs['kmeans5'] = kmeans.labels_.astype(str)
    '''
    kmeans = KMeans(n_clusters=15, random_state=0).fit(adata_spatial.X)
    adata_spatial.obs['kmeans15'] = kmeans.labels_.astype(str)
    kmeans = KMeans(n_clusters=5, random_state=0).fit(adata_spatial.X)
    adata_spatial.obs['kmeans5'] = kmeans.labels_.astype(str)
    kmeans = KMeans(n_clusters=10, random_state=0).fit(adata_spatial.X)
    adata_spatial.obs['kmeans10'] = kmeans.labels_.astype(str)
    '''
    print(adata_spatial)
    adata_spatial.write(modelname)
    return adata_spatial


