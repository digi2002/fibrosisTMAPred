#©2022 Qi Sun

import numpy as np
import matplotlib.pyplot as plt
from util import genadata
import pandas as pd
import scanpy as sc
from collections import Counter
from preprocessTMA import getOnsampleAdata,getOnsampleBatchCdata



def getOnsampeldata(modeldir,RES,datadir, fildid, samplefiles, modelname,keystr,labelDict,batch):
    if batch==False:
        adataonsample = getOnsampleAdata(datadir, fildid, samplefiles, modelname)
    elif batch==True:
        adataonsample=getOnsampleBatchCdata(modeldir, fildid, RES)
    pos = np.array(adataonsample.obsm['spatial'])
    xs=pos[:,0]
    ys=pos[:,1]
    df=adataonsample.to_df()
    nodefeatlabels=list(df.columns)
    nodefeats=df.values.tolist()
    xs_sel,ys_sel,nodefeatlabels,nodefeats_sel,labels_sel,\
    sexs_sel,histology_sel,stage_sel,age_sel,appalachia_sel=selDatabyCoor(xs, ys, nodefeatlabels, nodefeats, labelDict)

    if keystr=='sex':
        annotation_sel=sexs_sel
    elif keystr=='histology':
        annotation_sel=histology_sel
    elif keystr=='stage':
        annotation_sel=stage_sel
    elif keystr=='age':
        annotation_sel=age_sel
    elif keystr=='appalachia':
        annotation_sel = appalachia_sel
    elif keystr=='all':
        annotation_sel=[]
    return xs_sel, ys_sel, nodefeatlabels,nodefeats_sel, annotation_sel,labels_sel


def initialscatter(infile):
    num = 0
    from matplotlib.ticker import MultipleLocator
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1)
    xs, ys, nodefeatlabels, nodefeatlabelsdict, nodefeats = readOrifile(infile, num)
    ysnew=[]
    for y in ys:
        y=(-1)*y
        ysnew.append(y)
    ax1.scatter(xs, ysnew)
    plt.xticks(np.arange(-20.8, 3, 1.5))
    plt.yticks(np.arange(-8.6, 10, 1.5))
    user_interval = 1.5
    for _x in np.arange(-20.8, max(xs) + 1, user_interval):
        ax1.axvline(x=_x, ls='-', color='y')
    for _y in np.arange(-8.6, max(ys) + 1, user_interval):
        ax1.axhline(y=_y, ls='-')
    plt.show()
    print(len(xs))

def readOrifile(filename,num=0):
    nodefeatlabelsdict={}
    cnt=0
    f=open(filename)
    pre=f.readline()
    pre=f.readline()
    pre=f.readline()
    line = f.readline()
    tokens = line.strip().split('\t')[:25]
    nodefeatlabels = [float(i) for i in tokens]
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
            tokens=line.strip().split('\t')[:30]
            x=float(tokens[1])
            y=float(tokens[2])*(-1)
            feat=tokens[3:-2]
            nodefeat = [float(k) for k in feat]
            xs.append(x)
            ys.append(y)
            nodefeats.append(nodefeat)
    return np.array(xs),np.array(ys),nodefeatlabels,nodefeatlabelsdict,np.array(nodefeats)
'''
def selDatabyCoor(infile,labelDict):
    num=0
    xs, ys, nodefeatlabels, nodefeatlabelsdict, nodefeats = readOrifile(infile, num)
    xs_sel=[]
    ys_sel=[]
    labels_sel=[]
    nodefeats_sel=[]
    sexs_sel=[]
    histology_sel=[]
    stage_sel=[]
    age_sel=[]
    appalachia_sel=[]
    for x,y,nodefeat in zip(xs,ys,nodefeats):
        for label in labelDict.keys():
            x1,x2,y1,y2=labelDict[label]['x1'],labelDict[label]['x2'],labelDict[label]['y1'],labelDict[label]['y2']
            if x1<x<x2 and y1<y<y2:
                xs_sel.append(x)
                ys_sel.append(y)
                labels_sel.append(label)
                nodefeats_sel.append(nodefeat)
                sexs_sel.append(labelDict[label]['sex'])
                histology_sel.append(labelDict[label]['histology'])
                stage_sel.append(labelDict[label]['stage'])
                age_sel.append(labelDict[label]['age'])
                appalachia_sel.append(labelDict[label]['appalachia'])
                break
    print('# of original pixels: '+str(len(xs)))
    print('# of selected pixels: '+str(len(xs_sel)))
    nodefeats_sel=np.array(nodefeats_sel)
    return xs_sel,ys_sel,nodefeatlabels,nodefeatlabelsdict,nodefeats_sel,labels_sel,sexs_sel,histology_sel,stage_sel,age_sel,appalachia_sel
'''

def selDatabyCoor(xs, ys, nodefeatlabels, nodefeats,labelDict):
    for label in labelDict.keys():
        if 'race' in labelDict[label].keys():
            flag='vcu'
        else:
            flag='lung'
        break
    num=0
    xs_sel=[]
    ys_sel=[]
    labels_sel=[]
    nodefeats_sel=[]
    sexs_sel=[]
    histology_sel=[]
    stage_sel=[]
    age_sel=[]
    appalachia_sel=[]
    for x,y,nodefeat in zip(xs,ys,nodefeats):
        for label in labelDict.keys():
            x1,x2,y1,y2=labelDict[label]['x1'],labelDict[label]['x2'],labelDict[label]['y1'],labelDict[label]['y2']
            if x1<x<x2 and y1<y<y2:
                xs_sel.append(x)
                ys_sel.append(y)
                labels_sel.append(label)
                nodefeats_sel.append(nodefeat)
                sexs_sel.append(labelDict[label]['sex'])
                histology_sel.append(labelDict[label]['histology'])
                stage_sel.append(labelDict[label]['stage'])
                age_sel.append(labelDict[label]['age'])
                if flag=='vcu':
                    appalachia_sel.append(labelDict[label]['race'])
                else:
                    appalachia_sel.append(labelDict[label]['appalachia'])
                break
    print('# of on sample pixels: '+str(len(xs)))
    print('# of pixels with selected annotations: '+str(len(xs_sel)))
    nodefeats_sel=np.array(nodefeats_sel)
    return xs_sel,ys_sel,nodefeatlabels,nodefeats_sel,labels_sel,sexs_sel,histology_sel,stage_sel,age_sel,appalachia_sel




def labelscatter(infile,labelDict):
    num=0
    xs, ys, nodefeatlabels, nodefeatlabelsdict, nodefeats = readOrifile(infile, num)
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1)
    xs_sel=[]
    ys_sel=[]
    labels_sel=[]
    nodefeats_sel=[]
    sexs_sel=[]
    histology_sel=[]
    stage_sel=[]
    age_sel=[]
    appalachia_sel=[]
    for x,y,nodefeat in zip(xs,ys,nodefeats):
        for label in labelDict.keys():
            x1,x2,y1,y2=labelDict[label]['x1'],labelDict[label]['x2'],labelDict[label]['y1'],labelDict[label]['y2']
            if x1<x<x2 and y1<y<y2:
                xs_sel.append(x)
                ys_sel.append(y)
                labels_sel.append(label)
                nodefeats_sel.append(nodefeat)
                sexs_sel.append(labelDict[label]['sex'])
                histology_sel.append(labelDict[label]['histology'])
                stage_sel.append(labelDict[label]['stage'])
                age_sel.append(labelDict[label]['age'])
                appalachia_sel.append(labelDict[label]['appalachia'])
                break
    ax1.scatter(xs_sel, ys_sel,s=1)
    plt.xticks(np.arange(-20.8, 3, 1.5))
    plt.yticks(np.arange(-8.6, 10, 1.5))
    user_interval = 1.5
    for _x in np.arange(-20.8, max(xs) + 1, user_interval):
        ax1.axvline(x=_x, ls='-', color='y')
    for _y in np.arange(-8.6, max(ys) + 1, user_interval):
        ax1.axhline(y=_y, ls='-')
    plt.show()
    print('# of selected pixels: '+str(len(xs_sel)))
    nodefeats_sel=np.array(nodefeats_sel)
    return xs_sel,ys_sel,nodefeatlabels,nodefeatlabelsdict,nodefeats_sel,labels_sel,sexs_sel,histology_sel,stage_sel,age_sel,appalachia_sel

#Position	No.	Age	Sex	Organ/Anatomic Site	Pathology diagnosis	x1	x2	y1	y2

def checknodeanno(labelDict,keystr,keystrannos):
    newlabelDict={}
    for label in labelDict.keys():
        if labelDict[label][keystr] in keystrannos:
            newlabelDict[label]=labelDict[label]
    return newlabelDict

def readMultianno(annofile,patientnum):
    print('reading original the annotation file')
    if annofile.find('vcu')!=-1:
        flag='vcu'
    else:
        flag='lung'
    patientposes=[]
    sexes=[]
    histologies=[]
    stages=[]
    ages=[]
    appalachias=[]
    labelDict={}
    #f=open(annofile, mode='r',encoding='utf-8')
    f = open(annofile,encoding='utf-8')
    f.readline()
    for line in f:
        tokens=line.strip().split(',')
        if tokens[0]==patientnum:
            pos=tokens[1].strip()
            sex=tokens[2].strip()
            histology=tokens[3].strip()
            stage=tokens[4].strip()
            age=str(int(int(tokens[5].strip())/10))
            appalachia=tokens[6].strip()
            x1=float(tokens[7])
            x2=float(tokens[8])
            y1=float(tokens[9])
            y2=float(tokens[10])
            patientposes.append(patientnum+pos)
            sexes.append(sex)
            histologies.append(histology)
            stages.append(stage)
            ages.append(age)
            appalachias.append(appalachia)
            if flag=='vcu':
                labelDict[patientnum+pos]={'pos':patientnum+pos,'sex':sex,'histology':histology,
                                           'stage':stage,'age':age,'race':appalachia,'x1':x1,'x2':x2,'y1':y1,'y2':y2}
            elif flag == 'fibrosis':
                labelDict[patientnum + pos] = {'pos': pos, 'sex': sex, 'histology': histology,
                                               'type': stage, 'age': age, 'tissueid': appalachia, 'x1': x1, 'x2': x2,
                                               'y1': y1, 'y2': y2}
            else:
                labelDict[patientnum + pos] = {'pos': patientnum + pos, 'sex': sex, 'histology': histology,
                                               'stage': stage, 'age': age, 'appalachia': appalachia, 'x1': x1, 'x2': x2,
                                               'y1': y1, 'y2': y2}
            print('# of label: ' + str(len(patientposes)) + ' : ' + str(Counter(patientposes)))
            print('sex: ' + str(len(sexes)) + ' : ' + str(Counter(sexes)))
            print('histology: ' + str(len(histologies)) + ' : ' + str(Counter(histologies)))
            print('type: ' + str(len(stages)) + ' : ' + str(Counter(stages)))
            print('age: ' + str(len(ages)) + ' : ' + str(Counter(ages)))
            print('appalachia: ' + str(len(appalachias)) + ' : ' + str(Counter(appalachias)))
    return labelDict


def readAnno(annofile):
    print('reading original the annotation file')
    poses=[]
    sexes=[]
    histologies=[]
    stages=[]
    ages=[]
    appalachias=[]
    labelDict={}
    f=open(annofile, mode='r',encoding='utf-8-sig')
    f.readline()
    for line in f:
        tokens=line.strip().split(',')
        pos=tokens[0].strip()
        sex=tokens[1].strip()
        histology=tokens[2].strip()
        stage=tokens[3].strip()
        age=str(int(int(tokens[4].strip())/10))
        appalachia=tokens[5].strip()
        x1=float(tokens[6])
        x2=float(tokens[7])
        y1=float(tokens[8])
        y2=float(tokens[9])
        poses.append(pos)
        sexes.append(sex)
        histologies.append(histology)
        stages.append(stage)
        ages.append(age)
        appalachias.append(appalachia)
        labelDict[pos]={'pos':pos,'sex':sex,'histology':histology,'stage':stage,'age':age,'appalachia':appalachia,'x1':x1,'x2':x2,'y1':y1,'y2':y2}
    print('# of label: '+str(len(poses))+' : '+str(Counter(poses)))
    print('sex: '+ str(len(sexes))+ ' : '+str(Counter(sexes)))
    print('histology: '+str(len(histologies)) + ' : '+str(Counter(histologies)))
    print('stage: '+str(len(stages))+ ' : '+str(Counter(stages)))
    print('age: '+str(len(ages))+' : '+str(Counter(ages)))
    print('appalachia: '+str(len(appalachias))+' : '+str(Counter(appalachias)))
    return labelDict


def genadataAnno(labelDict,infile,fildid):
    xs_sel,ys_sel,nodefeatlabels,nodefeatlabelsdict,\
    nodefeats_sel,labels_sel,sexs_sel,histology_sel,stage_sel,age_sel,appalachia_sel=labelscatter(infile,labelDict)
    varindex = pd.DataFrame(index=nodefeatlabels)
    ys_selnew=[]
    for y in ys_sel:
        ys_selnew.append(y*(-1))
    adata = genadata(xs_sel, ys_selnew, nodefeats_sel, fildid, varindex)
    sc.pp.normalize_total(adata, target_sum=1, inplace=True)
    sc.pp.log1p(adata, base=2)
    adata.obs['label'] = labels_sel
    adata.obs['sex'] = sexs_sel
    adata.obs['histology'] = histology_sel
    adata.obs['stage'] = stage_sel
    adata.obs['age'] = age_sel
    adata.obs['appalachia'] = appalachia_sel
    print(adata)
    return adata


def genadataAnno(labelDict,adataonsample,fildid):
    for label in labelDict.keys():
        if 'race' in labelDict[label].keys():
            flag='vcu'
        else:
            flag='lung'
        break
    pos = np.array(adataonsample.obsm['spatial'])
    xs = pos[:, 0]
    ys = pos[:, 1]
    df = adataonsample.to_df()
    nodefeatlabels = list(df.columns)
    nodefeats = df.values.tolist()

    xs_sel, ys_sel, nodefeatlabels, nodefeats_sel, labels_sel, \
    sexs_sel, histology_sel, stage_sel, age_sel, appalachia_sel = selDatabyCoor(xs, ys, nodefeatlabels, nodefeats,
                                                                                labelDict)
    varindex = pd.DataFrame(index=nodefeatlabels)
    adata = genadata(xs_sel, ys_sel, nodefeats_sel, fildid, varindex)
    adata.obs['label'] = labels_sel
    adata.obs['sex'] = sexs_sel
    adata.obs['histology'] = histology_sel
    adata.obs['stage'] = stage_sel
    adata.obs['age'] = age_sel
    if flag=='vcu':
        adata.obs['race'] = appalachia_sel
    else:
        adata.obs['appalachia'] = appalachia_sel
    print(adata)
    return adata
