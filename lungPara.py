#©2022 Qi Sun


samplefiles = {'lung1B': '1B lung TMA 07082020new.txt',
               'lung2B': '2B lung TMA 07082020new.txt',
               'lung3B': '3B lung TMA 07092020new.txt',
               'lung4B': '4B lung TMA 07072020new.txt'}

keystrSet = {'sex': ['F', 'M'],
             'appalachia': ['non-Appalachia', 'Appalachia'],
             'histology': ['Squamous cell carcinoma', 'Adenocarcinoma'],
             'stage': ['IB', 'IIB', 'IA', 'IIIA'],
             'age': ['5', '6', '7']}

labelIDdict={'Squamous cell carcinoma':0,'Adenocarcinoma':1}
idLabeldict={0:'Squamous cell carcinoma',1:'Adenocarcinoma'}


cellDE={'tumor':['1743.5917','1581.5356','1419.4812'],
        'necrosis':['2028.7228','1136.4036'],
        'fibroblasts':['1809.6495','1444.5138','2057.7407','2174.7847']}
topsall=['1743.5917','1581.5356','1419.4812','2028.7228',
         '1136.4036','1809.6495','1444.5138','2057.7407','2174.7847']


RES=0.1

runloc = 'LOC'
if runloc == 'LOC':
    dir = './'
else:
    dir='/scratch/qsu226/TMAclassification/'
datadir = dir + 'data/lung/'
modeldir = dir + 'models_lung/'
figdir = dir + 'figures_lung/'
resdir=dir+'results_lung/'
metadatadir=datadir+'metadata/'
coredatadir=datadir+'coredata/'


annofile=datadir+'/multiannotation.csv'
matrixfile=dir+'data/matrix/matrix.txt'

def getAnnofile():
    return annofile

def getmatrixfile():
    return matrixfile

def getPara():
    return samplefiles,keystrSet,RES

def getDir():
    return datadir,modeldir,figdir,resdir,metadatadir,coredatadir

def getlabelIDdict():
    return labelIDdict

def getidLabeldict():
    return idLabeldict

def getCellDE():
    return cellDE,topsall