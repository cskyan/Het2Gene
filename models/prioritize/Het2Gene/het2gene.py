import argparse
import os
import pickle
import numpy

def outRank(out_dir,sorted_scores:dict,topn):
    # make dir
    os.makedirs(out_dir, exist_ok=True)
    # write result
    resultPath = os.path.join(out_dir, 'result.txt')
    with open(resultPath, 'w') as f:
        f.write("Rank")
        f.write("\t")
        f.write("EntrezID")
        f.write("\t")
        f.write("Symbol")
        f.write("\t")
        f.write("Score")
        f.write("\n")
        i = 1
        for k in sorted_scores:
            ENid = k[0]
            gSymbol = k[1]
            f.write(str(i))
            f.write("\t")
            f.write(ENid)
            f.write("\t")
            f.write(gSymbol)
            f.write("\t")
            f.write(str(sorted_scores[k]))
            f.write("\n")
            i+=1
            if i > topn :
                break


def get_key(d:dict,x) :
    for k in d :
        if d[k] == x:
            return k
    return "NotFound*"

path_emb_gene ="../../../data/embeddings/Het2Gene/embedding_gene.pickle"
path_emb_hpo = "../../../data/embeddings/Het2Gene/embedding_hpo.pickle"
path_map_gene = "../../../data/embeddings/Het2Gene/geneMap.pickle"
path_map_hpo = "../../../data/embeddings/Het2Gene/hpoMap.pickle"
path_AllGeneMap = "../../../data/embeddings/Het2Gene/AllGeneMap.pickle"  # symbol → ID
with open(path_emb_hpo, "rb+") as f:
    embedding_hpo = pickle.load(f)

with open(path_emb_gene, "rb+") as f:
    embedding_gene = pickle.load(f)

with open(path_map_hpo, "rb+") as f:
    hpo_map = dict(pickle.load(f))

with open(path_map_gene, "rb+") as f:
    geneMap = dict(pickle.load(f))

with open(path_AllGeneMap,"rb+") as f :
    AllGeneMap = dict(pickle.load(f))

parser = argparse.ArgumentParser(description='args for prioritize')
parser.add_argument('--hps', type=str, default="HP:0001511,HP:0002194,HP:0000750,HP:0000343,HP:0009748")
parser.add_argument('--out_dir',type=str,default="RankResult")
parser.add_argument('--topn',type=int,default=1000)
args=parser.parse_args()
hps = args.hps
out_dir = args.out_dir
topn = args.topn
hpSet = str(hps).split(",")

def prioritize(hpSet:list,topn:int,out_dir) :
    agg_hpo = numpy.zeros(300)
    for h in hpSet :
        try:
            #agg_hpo = agg_hpo + hpoWeight[str(h).replace(":", "_")] * numpy.array(embedding_hpo[int(hpo_map[h])])
            agg_hpo = agg_hpo +  numpy.array(embedding_hpo[int(hpo_map[h])])
        except :
            print("{} not found, automatically skipped".format(h))
            continue
        scores = {}
        for g in range(0, len(geneMap)):
            gemb = embedding_gene[g]
            score = numpy.dot((agg_hpo).tolist(), gemb)
            ENid = str(get_key(geneMap,g)).split("Entrez:")[-1]
            gSymbol = get_key(AllGeneMap, ENid)
            scores.update({(ENid,gSymbol):round(score,3)})
        sorted_scores = dict(sorted(scores.items(), key=lambda x: x[1], reverse=True))
        outRank(out_dir,sorted_scores,topn)

if __name__ == "__main__" :
    if len(hpSet) == 0 :
        print("There is an error in the HP set and input. Please check if the HP term is correct. If the input is correct, please confirm if the HP term is outdated")
    else:
        prioritize(hpSet, topn, out_dir)
        print("GeneRankList has been successfully generated!")
