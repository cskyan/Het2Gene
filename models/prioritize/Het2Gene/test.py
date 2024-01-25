import os
import pickle

import numpy

def get_ranking(x,ls:list):
    ls.sort()
    r = ls.index(x)
    return len(ls) - r


path_emb_gene ="../../../data/embeddings/wHet2Gene/embedding_gene.pickle"
path_emb_hpo = "../../../data/embeddings/wHet2Gene/embedding_hpo.pickle"
path_map_gene = "../../../data/embeddings/wHet2Gene/geneMap.pickle"
path_map_hpo = "../../../data/embeddings/wHet2Gene/hpoMap.pickle"
path_AllGeneMap = "../../../data/embeddings/wHet2Gene/AllGeneMap.pickle"  # symbol → ID
patientPath = "../../../data/patient/"

with open(path_emb_hpo, "rb+") as f:  # format :   gId,gIndex
    embedding_hpo = pickle.load(f)

with open(path_emb_gene, "rb+") as f:  # format :   gId,gIndex
    embedding_gene = pickle.load(f)

with open(path_map_hpo, "rb+") as f:
    hpo_map = dict(pickle.load(f))

with open(path_map_gene, "rb+") as f:  # 'Entrez:51663': 4211
    geneMap = dict(pickle.load(f))

with open(path_AllGeneMap,"rb+") as f :        # 'BHLHE22': '27319'
    AllGeneMap = dict(pickle.load(f))

with open("causalgeneList.pickle","rb+") as f :
    causalgeneList = list(pickle.load(f))



hitList = numpy.zeros(52)
MRR = 0
for fileName in os.listdir(patientPath) :
    cgName = causalgeneList[int(fileName.split("_")[-1])-1]
    hpSet = []

    agg = numpy.zeros(300)
    with open(os.path.join(patientPath,fileName),"r") as f :
        hp = (f.readline().replace("\n",""))
        #print(hp)
        while hp :
            try:
                agg = agg + numpy.array(embedding_hpo[int(hpo_map[hp])])
            except:
                print("{} not found, automatically skipped".format(hp))
                hp = f.readline().replace("\n", "")
                continue
            hp = f.readline().replace("\n","")

    scores = []
    agg_cg_dp = numpy.dot(agg.tolist(),embedding_gene[int(geneMap["Entrez:"+AllGeneMap[cgName]])])
    for g in range(0, len(geneMap)):
        gemb = embedding_gene[g]
        score = numpy.dot((agg).tolist(), gemb)
        scores.append(score)

    rank = get_ranking(agg_cg_dp,scores)
    MRR += (1 / rank)
    print(rank)
    if rank > 50 :
        continue
    else:
        for i in range(rank,51) :
            hitList[i] +=1


print("hit1:{}".format(hitList[1]/209))
print("hit5:{}".format(hitList[5]/209))
print("hit10:{}".format(hitList[10]/209))
print("hit20:{}".format(hitList[20]/209))
print("hit30:{}".format(hitList[30]/209))
print("hit40:{}".format(hitList[40]/209))
print("hit50:{}".format(hitList[50]/209))
print("MRR:{}".format(MRR/209))

