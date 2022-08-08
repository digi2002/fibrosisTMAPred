from util import readOrifile

samplefiles = {'lung1B': 'Fibrosis TMA pngiso June0122_1809_2.txt'}

keystrSet = {'stage': ['AT', 'NAT','Normal'],
             'sex':['M','F']}

labelIDdict={'WHITE':0,'BLACK OR AFRICAN AMERICAN':1}
idLabeldict={0:'WHITE',1:'BLACK OR AFRICAN AMERICAN'}


cellDE={'tumor':['1743.5917','1581.5356','1419.4812'],
        'necrosis':['2028.7228','1136.4036'],
        'fibroblasts':['1809.6495','1444.5138','2057.7407','2174.7847']}


topsall=['1954.6842','2101.7661','1911.6249','1663.5862','1976.6971',
         '1581.5356','1905.6450','1743.5917','1403.4844']



RES=0.1

runloc = 'LOC'
if runloc == 'LOC':
    dir = './'
else:
    dir='/scratch/qsu226/TMAclassification/'
datadir = dir + 'data/fibrosisTMA/'
modeldir = dir + 'models_fibrosisTMA/'
figdir = dir + 'figures_fibrosisTMA/'
resdir=dir+'results_fibrosisTMA/'
metadatadir=dir+'metadata_fibrosisTMA/'
coredatadir=datadir+'coredata/'


annofile=metadatadir+'fibrosistma.csv'
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


def getNodefeatlabels(fildid):
    infile = datadir + samplefiles[fildid]
    num = 1
    print('infile: ' + infile)
    xs, ys, nodefeatlabels, nodefeatlabelsdict, nodefeats = readOrifile(infile, num)
    return nodefeatlabels
