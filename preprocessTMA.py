import scanpy as sc
from util import readOrifile,genadata
import pandas as pd
from plot_util import plotsingle,plotsingleglycan
import matplotlib.pyplot as plt
import numpy as np
from tmaPara import getPara,getDir,getAnnofile
from anno_lungTMA import readMultianno,genadataAnno



def initialscatter(adataonsample):
    pos=np.array(adataonsample.obsm['spatial'])
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.scatter(pos[:,0], pos[:,1],s=1)
    xmin=min(pos[:,0])
    xmax=-0.2
    ymin=min(pos[:,1])
    ymax=9

    plt.xticks(np.arange(xmin,xmax, 2.5))
    plt.yticks(np.arange(ymin, ymax, 2.5))
    user_interval = 1
    for _x in np.arange(xmin, xmax + 1, user_interval):
        ax1.axvline(x=_x, ls='-', color='y')
    for _y in np.arange(ymin, ymax + 1, user_interval):
        ax1.axhline(y=_y, ls='-')

    plt.show()

def getOnsampleAdata(datadir,fildid,samplefiles,modelname):
    infile = datadir + samplefiles[fildid]
    num=0
    print('infile: '+infile)
    xs, ys, nodefeatlabels, nodefeatlabelsdict, nodefeats = readOrifile(infile, num)
    nodefeatlabels = [float(i) for i in nodefeatlabels]
    nodefeatlabels = [str(i) for i in nodefeatlabels]
    varindex = pd.DataFrame(index=nodefeatlabels)
    adata = genadata(xs, ys, nodefeats, fildid, varindex)
    sc.pp.normalize_total(adata, target_sum=1, inplace=True)
    sc.pp.log1p(adata, base=2)

    adatalabel = sc.read(modelname)
    #adata = sc.read(modelname)
    if datadir.find('lung')!=-1:
        adata.obs['leiden']=adatalabel.obs['leiden']
        adata.obsm = adatalabel.obsm
        labellist=adata.obs['leiden'].unique().tolist()
        labellist=[i for i in labellist if i!='0']
        adataonsample=adata[adata.obs['leiden'].isin(labellist), :]
    elif datadir.find('fibrosisTMA')!=-1:
        adata.obs['kmeans7'] = adatalabel.obs['kmeans7']
        adata.obsm = adatalabel.obsm
        labellist = adata.obs['kmeans7'].unique().tolist()
        labellist = [i for i in labellist if i != '0']
        adataonsample = adata[adata.obs['kmeans7'].isin(labellist), :]
    else:
        adata.obs['kmeans2'] = adatalabel.obs['kmeans2']
        adata.obsm = adatalabel.obsm
        labellist = adata.obs['kmeans2'].unique().tolist()
        labellist = [i for i in labellist if i != '1']
        adataonsample = adata[adata.obs['kmeans2'].isin(labellist), :]
    return adataonsample


def getTMAOnsampleAdata(datadir,modelname):
    adata = sc.read(modelname)
    if datadir.find('vcu')==-1:
        labellist=adata.obs['leiden'].unique().tolist()
        labellist=[i for i in labellist if i!='0']
        adataonsample=adata[adata.obs['leiden'].isin(labellist), :]
    elif datadir.find('fibrosisTMA')==1:
        labellist = adata.obs['kmeans7'].unique().tolist()
        labellist = [i for i in labellist if i != '0']
        adataonsample = adata[adata.obs['kmeans7'].isin(labellist), :]
    else:
        labellist = adata.obs['kmeans2'].unique().tolist()
        labellist = [i for i in labellist if i != '1']
        adataonsample = adata[adata.obs['kmeans2'].isin(labellist), :]
    return adataonsample



def getOnsampleBatchCdata(modeldir,fildid,RES):
    modelname = modeldir + 'alllung' + '_R' + str(RES) + '_onsample_batchcorrect.h5ad'
    adata = sc.read(modelname)
    adata = adata[adata.obs['library_id']==fildid, :]
    return adata

def getTMAOnsample(modeldir,datadir,samplefiles,annofile,fildid,RES):
    print('read orignial annotations')
    labelDict = readMultianno(annofile, fildid.replace('lung', ''))
    # 读入onsample pixels
    inmodel = modeldir + fildid + '_R' + str(RES) + '.h5ad'
    adataonsample = getOnsampleAdata(datadir, fildid, samplefiles, inmodel)
    # 对onsample pixel进行标注
    adata = genadataAnno(labelDict, adataonsample, fildid)
    return adata


if __name__ == "__main__":
    datadir, modeldir, figdir, resdir, metadatadir, coredatadir = getDir()
    samplefiles, keystrSet, RES = getPara()

#######################################################
    #adataonsample = getOnsampleAdata(datadir, fildid, samplefiles, modelname)
#########################################################
    fildid = 'lung1B'
    modelname = modeldir + fildid + '_R' + str(RES) + '.h5ad'
    adata = sc.read(modelname)
    # 打印数据坐标
    #initialscatter(adata)

    #
    #initialscatter(adata)
    #(adata, fildid, RES, figdir, 'leiden', 'allpixel')
    #plotsingle(adataonsample, fildid, RES, figdir, 'leiden', 'onsample')
    # 查看去除matrix后的数据

    #for glycan in ['2466.9033','2101.7661','2779.9678','1485.5420','2393.8560','1905.6450']: #'1485.542 ',
   # for glycan in ['1282.4606', '1136.4036', '1444.5138', '1298.4552', '2101.7661']:

    annofile = getAnnofile()
    ad = getTMAOnsample(modeldir, datadir, samplefiles, annofile, fildid, RES)
    #for glycan in ['1282.4606', '1136.4036', '1444.5138', '1298.4552', '2101.7661']:
    for glycan in ['1282.4606', '1136.4036', '1298.4552', '1485.542', '1444.5138']:

    #for glycan in ['1444.5138','1809.6495','1175.3761']:
        #for glycan in ['1809.6495', '1444.5138', '1791.6316', '2122.7349', '2057.7407', '1793.6561']:
        #plotsingle(adataonsample, fildid, RES, figdir, leiden, 'onsample')
        plotsingleglycan(ad, fildid,figdir, glycan,'rmmatrix')

        # ['1282.4606', '1136.4036', '1298.4552', '1485.542', '1444.5138']
        # all ['1282.4606', '1136.4036', '1298.4552', '1485.542', '1444.5138']

        # remove matrix['1282.4606', '1136.4036', '1444.5138', '1298.4552', '2101.7661']



