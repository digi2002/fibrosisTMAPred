#©2022 Qi Sun


import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np


def changelabel(adata,color):
    labelcountdict={}
    labelcountliststr=[]
    labellist = adata.obs[color].tolist()
    for i in labellist:
        if i not in labelcountdict.keys():
            labelcountdict[i]=0
        labelcountdict[i]=labelcountdict[i]+1
    labelcountlist=sorted(labelcountdict.items(), key=lambda x: x[1], reverse=True)
    for (label,count) in labelcountlist:
        labelcountliststr.append(label+':'+str(count))
    return labelcountliststr

def codeColorlist(colorlist):
    cnt=40
    newcolorlsit=[]
    for color in colorlist:
        newcolorlsit.append(str(cnt))
        cnt=cnt+1
    return newcolorlsit



def plotlabelglycan(adata,fildid,figdir,label,glycan,featcluster,app=''):
    print(adata)
    color=glycan
    markerrange = adata[:, adata.var_names == color].X.tolist()
    minv = np.quantile(markerrange, 0.1)
    maxv = np.quantile(markerrange, 0.9)
    fig, axs = plt.subplots(1, 2, figsize=(18, 7))
    sc.pl.spatial(
        adata,
        img_key=None,
        vmin=minv,
        vmax=maxv,
        library_id=None,
        color=str(color),
        title=str(color)+ featcluster[glycan],
        size=10,
        show=False,
        cmap = 'RdYlBu_r',
        ax=axs[1]

    )
    sc.pl.spatial(
        adata,
        img_key=None,
        library_id=None,
        color=label,
        title=label,
        size=10,
        show=False,
        cmap='RdYlBu_r',
        ax=axs[0]

    )
    plt.savefig(figdir+fildid+'_'+color+'_'+label+app+'.png')


def plotsingleglycan(adata,fildid,figdir,glycan,app=''):
    print(adata)
    color=glycan
    markerrange = adata[:, adata.var_names == color].X.tolist()
    #minv=0
    #maxv=1
    minv = np.quantile(markerrange, 0.01)
    maxv = np.quantile(markerrange, 0.99)
    sc.pl.spatial(
        adata,
        img_key=None,
        vmin=minv,
        vmax=maxv,
        library_id=None,
        color=str(color),
        title=str(color),
        size=10,
        show=False,
        cmap = 'RdYlBu_r',
    )
    plt.savefig(figdir+fildid+'_'+color+app+'.png')


def plotCorreanno(adata,type1,type2,glycan1,glycan2,glycan1type,glycan2type,figdir,fildid,RES,app):
    clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
    #colorlist=codeColorlist(adata.obs[label].unique().tolist())
    #fig, axs = plt.subplots(2, 2, figsize=(18, 7))
    fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(nrows=2, ncols=2,figsize=(18, 7))
    sc.pl.spatial(
        adata,
        img_key=None,
        # legend_loc='on data',
        # legend_fontsize=3,
        library_id=None,
        color=type1,
        title=type1,
        size=20,
        palette=[
            v
            for k, v in clusters_colors.items()
            if k in adata.obs[type1].unique().tolist()
        ],
        ax=ax1,
        show=False
    )

    markerrange = adata[:, adata.var_names == glycan1].X.tolist()
    minv = np.quantile(markerrange, 0.01)
    maxv = np.quantile(markerrange, 0.99)
    sc.pl.spatial(
        adata,
        img_key=None,
        #legend_loc='on data',
        #legend_fontsize=3,
        vmin=minv,
        vmax=maxv,
        library_id=None,
        color=glycan1,
        title=glycan1+' for '+glycan1type,
        size=20,
        ax=ax2,
        cmap='RdYlBu_r',
        show=False
    )

    sc.pl.spatial(
        adata,
        img_key=None,
        # legend_loc='on data',
        # legend_fontsize=3,
        library_id=None,
        color=type2,
        title='race',
        size=20,
        palette=[
            v
            for k, v in clusters_colors.items()
            if k in adata.obs[type2].unique().tolist()
        ],
        ax=ax3,
        show=False
    )
    #if label=='label':
    #    position='on data'
    #else:
    #    position='right margin'
    markerrange = adata[:, adata.var_names == glycan2].X.tolist()
    minv = np.quantile(markerrange, 0.01)
    maxv = np.quantile(markerrange, 0.99)
    sc.pl.spatial(
        adata,
        img_key=None,
        legend_loc='on data',
        legend_fontsize=5,
        library_id=None,
        vmin=minv,
        vmax=maxv,
        color=glycan2,
        title=glycan2+' for '+glycan2type,
        size=20,
        ax=ax4,
        cmap='RdYlBu_r',
        show=False
    )
    #print(supstr)
    #fig.suptitle(supstr)
    #plt.show()
    plt.savefig(figdir+fildid+'_'+glycan1+'_'+glycan2+'_R'+str(RES)+app+'.png')


def plotanno(adata,label,figdir,fildid,RES,app):
    clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
    colorlist=codeColorlist(adata.obs[label].unique().tolist())

    cluster='leiden'
    fig, axs = plt.subplots(1, 4, figsize=(18, 7))
    sc.pl.umap(
        adata, legend_loc='on data',color=[cluster], palette=[
            v
            for k, v in clusters_colors.items()
            if k in adata.obs[cluster].unique().tolist()],
        ax=axs[0],
        show=False,
    )

    labelcountlist=changelabel(adata,cluster)
    supstr=', '.join(labelcountlist)
    sc.pl.spatial(
        adata,
        img_key=None,
        #legend_loc='on data',
        #legend_fontsize=3,
        library_id=None,
        color=cluster,
        title=fildid,
        size=20,
        palette=[
            v
            for k, v in clusters_colors.items()
            if k in adata.obs[cluster].unique().tolist()
        ],
        ax=axs[1],
        show=False
    )

    sc.pl.umap(
        adata, legend_loc='on data', color=[label], palette=[
            v
            for k, v in clusters_colors.items()
            if k in colorlist],
        ax=axs[2],
        show=False,
    )

    labelcountlist = changelabel(adata, label)
    supstr = ', '.join(labelcountlist)

    if label=='label':
        position='on data'
    else:
        position='right margin'
    sc.pl.spatial(
        adata,
        img_key=None,
        legend_loc=position,
        legend_fontsize=5,
        library_id=None,
        color=label,
        title=fildid,
        size=20,
        palette=[
            v
            for k, v in clusters_colors.items()
            if k in colorlist
        ],
        ax=axs[3],
        show=False
    )

    print(supstr)
    fig.suptitle(supstr)
    #plt.show()
    plt.savefig(figdir+fildid+'_'+label+'_R'+str(RES)+app+cluster+'.png')



def plotannotmp(adata,label,figdir,fildid,RES,app):
    clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
    colorlist=codeColorlist(adata.obs[label].unique().tolist())

    #fig, axs = plt.subplots(1, 2, figsize=(18, 7))

    if label=='label':
        position='on data'
    else:
        position='right margin'
    sc.pl.spatial(
        adata,
        img_key=None,
        legend_loc=position,
        legend_fontsize=5,
        library_id=None,
        color=label,
        title=fildid,
        size=20,
        palette=[
            v
            for k, v in clusters_colors.items()
            if k in colorlist
        ],
        #ax=axs[1],
        show=False
    )
    plt.savefig(figdir+fildid+'_'+label+'_R'+str(RES)+app+'.png')


def classifyResults(siglabel,adata_train,adata_test_answer,adata_test_result,adata_test_coreresult,annotation_sel,correct,wrong,acc,figdir,keystr,classifier):
    clusters_colors = dict(zip([i for i in sorted(set(annotation_sel))], sc.pl.palettes.default_102))
    fig, axs = plt.subplots(1, 5, figsize=(18, 7))


    color = 'label'
    sc.pl.spatial(
        adata_train,
        img_key=None,
        # legend_loc='on data',
        # legend_fontsize=3,
        library_id=None,
        color=color,
        title='train data',
        size=20,
        palette=[
            v
            for k, v in clusters_colors.items()
            if k in codeColorlist(adata_train.obs['label'].unique().tolist())
        ],
        ax=axs[0],
        show=False
    )
    title = 'No1 DE: ' + str(siglabel)
    feat = str(siglabel)
    markerrange = adata_train[:, adata_train.var_names == feat].X.tolist()
    minv = np.quantile(markerrange, 0.01)
    maxv = np.quantile(markerrange, 0.99)

    sc.pl.spatial(
        adata_train,
        img_key=None,
        # legend_loc='on data',
        # legend_fontsize=3,
        library_id=None,
        color=feat,
        title=title,
        size=20,
        palette=[
            v
            for k, v in clusters_colors.items()
            if k in codeColorlist(adata_train.obs['label'].unique().tolist())
        ],
        ax=axs[1],
        show=False
    )

    sc.pl.spatial(
        adata_test_answer,
        img_key=None,
        # legend_loc='on data',
        # legend_fontsize=3,
        library_id=None,
        color=color,
        title='ground truth',
        size=20,
        palette=[
            v
            for k, v in clusters_colors.items()
            if k in codeColorlist(adata_test_answer.obs['label'].unique().tolist())
        ],
        ax=axs[2],
        show=False
    )

    title = 'pixel prediction: ' + str(acc)
    sc.pl.spatial(
        adata_test_result,
        img_key=None,
        # legend_loc='on data',
        # legend_fontsize=3,
        library_id=None,
        color=color,
        title=title,
        size=20,
        palette=[
            v
            for k, v in clusters_colors.items()
            if k in codeColorlist(adata_test_result.obs['label'].unique().tolist())
        ],
        ax=axs[3],
        show=False
    )

    coretitle = 'core prediction: ' + str(round(correct / (correct + wrong), 2)) + '\n correct: ' + str(
        correct) + ' wrong: ' + str(wrong)
    sc.pl.spatial(
        adata_test_coreresult,
        img_key=None,
        # legend_loc='on data',
        # legend_fontsize=3,
        library_id=None,
        color=color,
        title=coretitle,
        size=20,
        palette=[
            v
            for k, v in clusters_colors.items()
            if k in codeColorlist(adata_test_coreresult.obs['label'].unique().tolist())
        ],
        ax=axs[4],
        show=False
    )

    plt.savefig(figdir + keystr + '_core' + color + classifier+'.png')

def plotsingle(adata,sampleid,RES,figdir,color,app=''):
    print(sampleid)
    print(adata)
    clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))

    colorlist=adata.obs[color].unique().tolist()

    fig, axs = plt.subplots(1, 2, figsize=(12, 7))
    sc.pl.umap(
        adata, legend_loc='on data',color=[color], palette=[
            v
            for k, v in clusters_colors.items()
            if k in colorlist],
        ax=axs[0],
        show=False,
    )

    labelcountlist=changelabel(adata,color)
    supstr=', '.join(labelcountlist)
    sc.pl.spatial(
        adata,
        img_key=None,
        #legend_loc='on data',
        #legend_fontsize=3,
        library_id=None,
        color=color,
        title=sampleid,
        size=20,
        palette=[
            v
            for k, v in clusters_colors.items()
            if k in colorlist
        ],
        ax=axs[1],
        show=False
    )
    print(supstr)
    fig.suptitle(supstr)
    #plt.show()
    plt.savefig(figdir+sampleid+'_'+color+'_R'+str(RES)+app+'.png')



def umapplot(adata,color):
    clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
    sc.pl.umap(
        adata, color=[color,'library_id'], palette=[
            v
            for k, v in clusters_colors.items()
            if k in adata.obs.leiden.unique().tolist()],
    ) #legend_loc='on data',

def singlespatial(adata,sampleid,RES,figdir,color,title,app=''):
    print(sampleid)
    print(adata)
    #clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
    #colorlist=adata.obs[color].unique().tolist()

    clusters_colors = {'0_normal': 'tab:blue', '1_fibrosis': 'tab:orange', '2_nonfib': 'tab:green'}
    palette = [
        v
        for k, v in clusters_colors.items()
        if k in adata.obs[color].unique().tolist()]


    labelcountlist=changelabel(adata,color)
    supstr=', '.join(labelcountlist)
    sc.pl.spatial(
        adata,
        img_key=None,
        #legend_loc='on data',
        #legend_fontsize=3,
        library_id=None,
        color=color,
        title=title,
        size=20,
        palette = palette,
        #palette=[
        #    v
        #    for k, v in clusters_colors.items()
        #    if k in colorlist
        #],
        show=False
    )
    print(supstr)
    #plt.show()
    plt.savefig(figdir+sampleid+'_'+color+'_R'+str(RES)+app+'.png')


def plotsingleglycangrey(modeldir,sampleid,figdir,colorlist):
    color='leiden'
    adata = sc.read(modeldir+sampleid+'.h5ad')
    print(sampleid)
    print(adata)

    allcolor=adata.obs[color].unique().tolist()

    greycolorlist=[]
    for i in allcolor:
        if i not in colorlist:
            greycolorlist.append(i)

    indicegrey = []
    indice = []

    for feat in greycolorlist:
        ad = adata[adata.obs[color] == feat, :]
        indicegrey.extend(ad.obs.index)
    adgrey=adata[indicegrey, :]
    adgrey.obs[color]='NA'

    for feat in colorlist:
        ad = adata[adata.obs[color] == feat, :]
        indice.extend(ad.obs.index)
    ad=adata[indice, :]

    adata=adgrey.concatenate(ad)

    clusters_colors = dict(zip([str(i) for i in range(102)], sc.pl.palettes.default_102))
    clusters_colors_grey={}
    for cluster,colorname in clusters_colors.items():
        if cluster in colorlist:
            clusters_colors_grey[cluster] = colorname
    clusters_colors_grey['NA'] = 'lightgrey'

    fig, axs = plt.subplots(1, 2, figsize=(12, 7))
    sc.pl.umap(
        adata, legend_loc='on data',color=[color], palette=[
            v
            for k, v in clusters_colors_grey.items()
            if k in adata.obs[color].unique().tolist()],
        ax=axs[0],
        show=False
    )

    labelcountlist=changelabel(adata,color)
    supstr=', '.join(labelcountlist)
    sc.pl.spatial(
        adata,
        img_key=None,
        #legend_loc='on data',
        #legend_fontsize=3,
        library_id=None,
        color=color,
        title=sampleid,
        size=20,
        palette=[
            v
            for k, v in clusters_colors_grey.items()
            if k in adata.obs[color].unique().tolist()
        ],
        ax=axs[1],
        show=False
    )
    print(supstr)
    fig.suptitle(supstr)
    plt.savefig(figdir + sampleid + '_' + color + 'grey_cluster167.png')

