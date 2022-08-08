import pandas as pd
import scanpy as sc
from plot_util import plotsingleglycan
import pandas as pd

def readmatrix():
    matrixlist=[]
    f=open('./data/matrix/matrix.txt')
    for line in f:
        matrixlist.append(line.strip())
    return matrixlist

def DEana_foldchange(adata,fildid,clusterid,topN,featlist,all=False):
    matrixlist=readmatrix()
    if all==False:
        adata = adata[adata.obs['leiden'].isin(featlist), :]
    METHOD = "wilcoxon"
    sc.tl.rank_genes_groups(adata, method=METHOD, tie_correct=True,groups=[clusterid], groupby='anno', key_added='de')
    result = adata.uns['de']
    groups = [clusterid]
    sta = pd.DataFrame(
        {group + '_' + key: result[key][group]
         for group in groups for key in ['names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges']})
    name=clusterid + '_names'
    foldchange=clusterid + '_scores'#zhuyi shi zhege!!!
    inverse_boolean_series = ~sta[name].isin(matrixlist)
    sta = sta[inverse_boolean_series]
    sta=sta.sort_values(by=foldchange, ascending=False).head(topN)
    tops=sta[name].tolist()
    sta.to_csv('./heatmaptmp/' + fildid +'-'+clusterid+ '-top'+str(topN)+'.csv', sep=',')
    return tops



def deTissueROC(fildid,clusterid,topk):
    adata = sc.read('../covid/' + 'models_resolution1/' + fildid + 'DE.h5ad')
    print('../covid/' + 'models_resolution1/' + fildid + 'DE.h5ad')
    tops=DEana_foldchange(adata, fildid, clusterid, topk, [], all=True)
    #adata = adata[adata.obs['leiden']==clusterid, :]
    return tops,adata

def genTMAfeat(topN,ad):
    #ad = sc.read('./models_fibrosisTMA/lung1B_R0.1.h5ad')
    fildid = 'lung1B'

    normallabel = ['1BG1', '1BG2', '1BG3', '1BG4', '1BG5', '1BG6', '1BG7', '1BG8']
    fibrosislabel = ['1BF1', '1BF2', '1BF3', '1BF4', '1BF5', '1BF6', '1BF7', '1BF8',
                     '1BE1', '1BE2', '1BE3', '1BE4', '1BE5', '1BE6', '1BE7', '1BE8',
                     '1BD1', '1BD2', '1BD3', '1BD4', '1BD5', '1BD6', '1BD7', '1BD8',
                     '1BC1', '1BC2', '1BC3', '1BC4', '1BC5', '1BC6', '1BC7', '1BC8',
                     '1BB1', '1BB2', '1BB3', '1BB4', '1BB5', '1BB6', '1BB7', '1BB8',
                     '1BA1', '1BA2', '1BA3', '1BA4', '1BA5', '1BA6', '1BA7', '1BA8']

    adatanormal = ad[ad.obs['label'].isin(normallabel), :]
    adatafibrosis = ad[ad.obs['label'].isin(fibrosislabel), :]
    adatanormal.obs['anno'] = 'normal'
    adatafibrosis.obs['anno'] = 'fibrosis'
    adata = adatanormal.concatenate(adatafibrosis)

    tops = DEana_foldchange(adata, fildid, 'fibrosis', topN, [], all=True)
    return tops,adata


if __name__ == "__main__":
    fildid ='lung1B'
    topN =5
    tops, adata = genTMAfeat(topN)
    for glycan in tops:
        plotsingleglycan(adata, fildid, './figures_fibrosisTMA/', glycan, 'fibrosis')